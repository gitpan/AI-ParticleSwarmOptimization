package AI::ParticleSwarmOptimization;

use strict;
use warnings;
use Math::Random::MT qw(srand rand);

require Exporter;

our @ISA     = qw(Exporter);
our @EXPORT  = qw();
our $VERSION = '1.003';


sub new {
    my ($class, %params) = @_;
    my $self = bless {}, $class;

    $self->setParams (%params);
    srand ($self->{randSeed});
    return $self;
}


sub setParams {
    my ($self, %params) = @_;

    if (defined $params{-fitFunc}) {
        # Process required parameters - -fitFunc and -dimensions
        if ('ARRAY' eq ref $params{-fitFunc}) {
            ($self->{fitFunc}, @{$self->{fitParams}}) = @{$params{-fitFunc}}
        } else {
            $self->{fitFunc} = $params{-fitFunc};
        }

        $self->{fitParams} ||= [];
    }

    $self->{prtcls} = [] # Need to reinit if num dimensions changed
        if  defined $params{-dimensions}
        and defined $self->{dimensions}
        and $params{-dimensions} != $self->{dimensions};

    srand ($self->{randSeed}) # Need to reseed if seed changed
        if  defined $params{-randSeed}
        and defined $self->{randSeed}
        and $params{-randSeed} != $self->{randSeed};

    $self->{$_} = $params{"-$_"} for grep {exists $params{"-$_"}}
        qw/
        dimensions
        exitFit
        inertia
        iterations
        meWeight
        numNeighbors
        numParticles
        posMax
        posMin
        randSeed
        randStartVelocity
        stallSpeed
        themWeight
        verbose
        /;

    die "-dimensions must be greater than 0\n"
        unless ! exists $params{-dimensions} or $params{-dimensions} > 0;

    $self->{numParticles}   ||= $self->{dimensions} * 10 if defined $self->{dimensions};
    $self->{numNeighbors}   ||= int sqrt $self->{numParticles} if defined $self->{numParticles};
    $self->{iterations}     ||= 1000;
    $self->{posMax}         = 100 unless defined $self->{posMax};
    $self->{posMin}         = -$self->{posMax} unless defined $self->{posMin};
    $self->{meWeight}       ||= 0.5;
    $self->{themWeight}     ||= 0.5;
    $self->{inertia}        ||= 0.9;
    $self->{randSeed}       = rand (0xffffffff) unless defined $self->{randSeed};
    $self->{verbose}        ||= 0;

    return 1;
}


sub init {
    my ($self) = @_;

    die "-fitFunc must be set before init or optimize is called"
        unless $self->{fitFunc} and 'CODE' eq ref $self->{fitFunc};
    die "-dimensions must be set to 1 or greater before init or optimize is called"
        unless $self->{dimensions} and $self->{dimensions} >= 1;

    $self->{prtcls} = [];
    $self->{bestBest} = undef;
    $self->_initParticles ();
    $self->{iterCount} = 0;

    # Normalise weights.
    my $totalWeight = $self->{inertia} + $self->{themWeight} + $self->{meWeight};

    $self->{inertia}    /= $totalWeight;
    $self->{meWeight}   /= $totalWeight;
    $self->{themWeight} /= $totalWeight;

    die "-posMax must be greater than -posMin"
        unless $self->{posMax} > $self->{posMin};
    $self->{$_} > 0 or die "-$_ must be greater then 0"
        for qw/numParticles/;

    $self->{deltaMax} = ($self->{posMax} - $self->{posMin}) / 100.0;

    return 1;
}


sub optimize {
    my ($self, $iterations) = @_;

    $iterations ||= $self->{iterations};
    $self->init () unless $self->{prtcls};
    return $self->_swarm ($iterations);
}


sub getBestParticles {
    my ($self, $num) = @_;
    my @bests  = 0 .. $self->{numParticles} - 1;
    my $prtcls = $self->{prtcls};

    @bests = sort {$prtcls->[$a]{bestFit} <=> $prtcls->[$b]{bestFit}} @bests;
    $num ||= 1;
    return @bests[0 .. $num - 1];
}


sub getParticleBestPos {
    my ($self, $prtcl) = @_;

    return undef if $prtcl >= $self->{numParticles};
    $prtcl = $self->{prtcls}[$prtcl];

    return ($prtcl->{bestFit}, @{$prtcl->{bestPos}});
}


sub getIterationCount {
    my ($self) = @_;

    return $self->{iterCount};
}


sub _initParticles {
    my ($self) = @_;

    for my $id (0 .. $self->{numParticles} - 1) {
        $self->{prtcls}[$id]{id} = $id;
        $self->_initParticle ($self->{prtcls}[$id]);
    }
}


sub _initParticle {
    my ($self, $prtcl) = @_;

    # each particle is a hash of arrays with the array sizes being the
    # dimensionality of the problem space
    for my $d (0 .. $self->{dimensions} - 1) {
        $prtcl->{currPos}[$d] = _randInRange ($self->{posMin}, $self->{posMax});

        $prtcl->{velocity}[$d] =
            $self->{randStartVelocity}
            ? _randInRange (-$self->{deltaMax}, $self->{deltaMax})
            : 0;
    }

    $prtcl->{currFit} = $self->_calcPosFit ($prtcl->{currPos});
    $self->_calcNextPos ($prtcl);

    unless (defined $prtcl->{bestFit}) {
        $prtcl->{bestPos}[$_] = _randInRange ($self->{posMin}, $self->{posMax})
            for 0 .. $self->{dimensions} - 1;
        $prtcl->{bestFit} = $self->_calcPosFit ($prtcl->{bestPos});
    }
}


sub _calcPosFit {
    my ($self, $pos) = @_;

    return $self->{fitFunc}->(@{$self->{fitParams}}, @$pos);
}


sub _swarm {
    my ($self, $iterations) = @_;

    for my $iter (1 .. $iterations) {
        ++$self->{iterCount};
        last if defined $self->_moveParticles ($iter);

        $self->_updateVelocities ($iter);
    }

    return $self->{bestBest};
}


sub _moveParticles {
    my ($self, $iter) = @_;

    print "Iter $iter\n" if $self->{verbose} >= 2;

    for my $prtcl (@{$self->{prtcls}}) {
        @{$prtcl->{currPos}} = @{$prtcl->{nextPos}};
        $prtcl->{currFit} = $prtcl->{nextFit};

        my $fit = $prtcl->{currFit};

        if ($self->_betterFit ($fit, $prtcl->{bestFit})) {
            # Save position - best fit for this particle so far
            $self->_saveBest ($prtcl, $fit, $iter);
        }

        return $fit if defined $self->{exitFit} and $fit < $self->{exitFit};

        next unless $self->{verbose} >= 2;
        printf "Part %3d fit %8.2f", $prtcl->{id}, $fit
            if $self->{verbose} >= 2;
        printf " (%s @ %s)",
            join (', ', map {sprintf '%5.3f', $_} @{$prtcl->{velocity}}),
            join (', ', map {sprintf '%5.2f', $_} @{$prtcl->{currPos}})
            if $self->{verbose} >= 3;
        print "\n";
    }

    return undef;
}


sub _saveBest {
    my ($self, $prtcl, $fit, $iter) = @_;

    # for each dimension, set the best position as the current position
    @{$prtcl->{bestPos}} = @{$prtcl->{currPos}};

    $prtcl->{bestFit} = $fit;
    if ($self->_betterFit ($fit, $self->{bestBest})) {
        if ($self->{verbose} >= 1) {
            my $velSq;

            $velSq += $_ ** 2 for @{$prtcl->{velocity}};

            printf "#%05d: Particle $prtcl->{id} best: %.4f (vel: %.3f)\n",
                $iter, $fit, sqrt ($velSq)
        }

        $self->{bestBest} = $fit;
    }
}


sub _betterFit {
    my ($self, $new, $old) = @_;

    return ! defined ($old) || ($new < $old);
}


sub _updateVelocities {
    my ($self, $iter) = @_;

    for my $prtcl (@{$self->{prtcls}}) {
        my $bestN = $self->{prtcls}[$self->_getBestNeighbour ($prtcl)];
        my $velSq;

        for my $d (0 .. $self->{dimensions} - 1) {
            my $meFactor = _randInRange (-$self->{meWeight}, $self->{meWeight});
            my $themFactor =
                _randInRange (-$self->{themWeight}, $self->{themWeight});
            my $meDelta = $prtcl->{bestPos}[$d] - $prtcl->{currPos}[$d];
            my $themDelta = $bestN->{bestPos}[$d] - $prtcl->{currPos}[$d];

            $prtcl->{velocity}[$d] =
                  $prtcl->{velocity}[$d] * $self->{inertia}
                + $meFactor * $meDelta
                + $themFactor * $themDelta;
            $velSq += $prtcl->{velocity}[$d] ** 2;
        }

        if (! $velSq or $self->{stallSpeed} and $velSq <= $self->{stallSpeed}) {
            $self->_initParticle ($prtcl);
            printf "#%05d: Particle $prtcl->{id} stalled\n", $iter
                if $self->{verbose} > 1;
        }

        $self->_calcNextPos ($prtcl);
    }
}


sub _calcNextPos {
    my ($self, $prtcl) = @_;

    for my $d (0 .. $self->{dimensions} - 1) {
        $prtcl->{nextPos}[$d] = $prtcl->{currPos}[$d] + $prtcl->{velocity}[$d];
        if ($prtcl->{nextPos}[$d] < $self->{posMin}) {
            $prtcl->{nextPos}[$d]  = $self->{posMin};
            $prtcl->{velocity}[$d] = 0;
        } elsif ($prtcl->{nextPos}[$d] > $self->{posMax}) {
            $prtcl->{nextPos}[$d]  = $self->{posMax};
            $prtcl->{velocity}[$d] = 0;
        }
    }

    $prtcl->{nextFit} = $self->_calcPosFit ($prtcl->{nextPos});
}


sub _randInRange {
    my ($min, $max) = @_;
    return $min + rand ($max - $min);
}


sub _getBestNeighbour {
    my ($self, $prtcl) = @_;
    my $bestNFitness;
    my $bestNIndex;

    for my $neighbor (0 .. $self->{numNeighbors} - 1) {
        my $prtclNIndex = ($prtcl + $neighbor) % $self->{numParticles};

        if (!defined ($bestNFitness)
            || $self->{prtcls}[$prtclNIndex]{bestFit} < $bestNFitness)
        {
            $bestNFitness = $self->{prtcls}[$prtclNIndex]{bestFit};
            $bestNIndex   = $prtclNIndex;
        }
    }

    return $bestNIndex;
}


1;


=head1 NAME

AI::ParticleSwarmOptimization - Particle Swarm Optimization (object oriented)

=head1 SYNOPSIS

    use AI::ParticleSwarmOptimization;

    my $pso = AI::ParticleSwarmOptimization->new (
        fitFunc    => \&calcFit,
        dimensions => 3,
        );
    my $fitValue       = $pso->optimize ();
    my ($best)         = $pso->getBestParticles (1);
    my ($fit, @values) = $pso->getParticleBestPos ($best);

    printf "Fit %.4f at (%s)\n",
        $fit, join ', ', map {sprintf '%.4f', $_} @values;


    sub calcFit {
        my @values = @_;
        my $offset = int (-@values / 2);
        my $sum;

        $sum += ($_ - $offset++) ** 2 for @values;
        return $sum;
    }

=head1 Description

The Particle Swarm Optimization technique uses communication of the current best
position found between a number of particles moving over a hyper surface as a
technique for locating the best location on the surface (where 'best' is the
minimum of some fitness function).

This pure Perl module is an implementation of the Particle Swarm Optimization
technique for finding minima of hyper surfaces. It presents an object oriented
interface that facilitates easy configuration of the optimization parameters and
(in principle) allows the creation of derived classes to reimplement all aspects
of the optimization engine (a future version will describe the replaceable
engine components).

This implementation allows communication of a local best point between a
selected number of neighbours. It does not support a single global best position
that is known to all particles in the swarm.

=head1 Methods

AI::ParticleSwarmOptimization provides the following public methods. The parameter lists shown
for the methods denote optional parameters by showing them in [].

=over 4

=item new (%parameters)

Create an optimization object. The following parameters may be used:

=over 4

=item I<-dimensions>: positive number, required

The number of dimensions of the hypersurface being searched.

=item I<-exitFit>: number, optional

If provided I<-exitFit> specifies allows an early termination of optimize if the
fitness value becomes equal or less than I<-exitFit>.

=item I<-fitFunc>: required

I<-fitFunc> is a reference to the fitness function used by the search. If extra
parameters need to be passed to the fitness function and array ref may be used
with the code ref as the first array element and parameters to be passed into
the fitness function as following elements. User provided parameters are passed
as the first parameters to the fitness function when it is called:

    my $pso = AI::ParticleSwarmOptimization->new (
        fitFunc    => [\&calcFit, $context],
        dimensions => 3,
        );

    ...

    sub calcFit {
        my ($context, @values) = @_;
        ...
        return $fitness;
        }

In addition to any user provided parameters the list of values representing the
current particle position in the hyperspace is passed in. There is one value per
hyperspace dimension.

=item I<-inertia>: positive or zero number, optional

Determines what proportion of the previous velocity is carried forward to the
next iteration. Defaults to 0.9

See also I<-meWeight> and I<-themWeight>.

=item I<-iterations>: number, optional

Number of optimization iterations to perform. Defaults to 1000.

=item I<-meWeight>: number, optional

Coefficient determining the influence of the current local best position on the
next iterations velocity. Defaults to 0.5.

See also I<-inertia> and I<-themWeight>.

=item I<-numNeighbors>: positive number, optional

Number of local particles considered to be part of the neighbourhood of the
current particle. Defaults to the square root of the total number of particles.

=item I<-numParticles>: positive number, optional

Number of particles in the swarm. Defaults to 10 times the number of dimensions.

=item I<-posMax>: number, optional

Maximum coordinate value for any dimension in the hyper space. Defaults to 100.

=item I<-posMin>: number, optional

Minimum coordinate value for any dimension in the hyper space. Defaults to
-I<-posMax> (if I<-posMax> is negative I<-posMin> should be set more negative).

=item I<-randSeed>: number, optional

Seed for the random number generator. Useful if you want to rerun an
optimization, perhaps for benchmarking or test purposes.

=item I<-randStartVelocity>: boolean, optional

Set true to initialize particles with a random velocity. Otherwise particle
velocity is set to 0 on initalization.

A range based on 1/100th of -I<-posMax> - I<-posMin> is used for the initial
speed in each dimension of the velocity vector if a random start velocity is
used.

=item I<-stallSpeed>: positive number, optional

Speed below which a particle is considered to be stalled and is repositioned to
a new random location with a new initial speed.

By default I<-stallSpeed> is undefined but particles with a speed of 0 will be
repositioned.

=item I<-themWeight>: number, optional

Coefficient determining the influence of the neighbourhod best position on the
next iterations velocity. Defaults to 0.5.

See also I<-inertia> and I<-meWeight>.

=item I<-verbose>: number, optional

If set to a non-zero value I<-verbose> determines the level of diagnostic print
reporting that is generated during optimization.

=back

=item B<setParams (%parameters)>

Set or change optimization parameters. See I<-new> above for a description of the
parameters that may be supplied.

=item B<init ()>

Reinitialize the optimization. B<init ()> will be called during the first call
to B<optimize ()> if it hasn't already been called.

=item B<optimize ()>

Runs the minimization optimization. Returns the fit value of the best fit
found. The best possible fit is negative infinity.

=item B<getParticleState ()>

Returns the vector of position

=item B<getBestParticles ([$n])>

Takes an optional count.

Returns a list containing the best $n prtcl numbers. If $n is not specified only
the best prtcl number is returned.

=item B<getParticleBestPos ($particleNum)>

Returns a list containing the best value of the fit and the vector of its point
in hyper space.

    my ($fit, @vector) = $pso->getParticleBestPos (3)

=item B<getIterationCount ()>

Return the number of iterations performed. This may be useful when the
I<-exitFit> criteria has been met or where multiple calls to I<optimize> have
been made.

=back

=head1 BUGS

Please report any bugs or feature requests to
C<bug-AI-PSO-OO at rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=AI-PSO-OO>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 SUPPORT

This module is supported by the author through CPAN. The following links may be
of assistance:

=over 4

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/AI-PSO-OO>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/AI-PSO-OO>

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=AI-PSO-OO>

=item * Search CPAN

L<http://search.cpan.org/dist/AI-PSO-OO>

=back

=head1 SEE ALSO

http://en.wikipedia.org/wiki/Particle_swarm_optimization

=head1 ACKNOWLEDGEMENTS

This module is an evolution of the AI::PSO module created by Kyle Schlansker.

=head1 AUTHOR

    Peter Jaquiery
    CPAN ID: GRANDPA
    grandpa@cpan.org

=head1 COPYRIGHT AND LICENSE

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=cut
