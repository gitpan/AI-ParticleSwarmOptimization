use strict;
use warnings;
use Test::More;

=head1 NAME

AI::ParticleSwarmOptimization test suite

=head1 DESCRIPTION

Test AI::ParticleSwarmOptimization

=cut

BEGIN {
    use lib '../lib'; # For development testing
    use AI::ParticleSwarmOptimization;

    plan (tests => 27);

    use_ok ("AI::ParticleSwarmOptimization");
}

ok (my $pso = AI::ParticleSwarmOptimization->new (), 'Constructor');
mustDie ('$pso->setParams (-fitFunc => 1)', 'Bad -fitFunc');
ok ($pso->setParams (-fitFunc => \&fitFunc,), 'Good -fitFunc (setParams)');
ok ($pso = AI::ParticleSwarmOptimization->new (-fitFunc => \&fitFunc,), 'Good -fitFunc (new)');
ok ($pso->setParams (-fitFunc => [\&fitFunc, 1]), 'Good -fitFunc (array)');

mustDie ('$pso->setParams (-dimensions => 0)', '-dimensions 0');
mustDie ('$pso->setParams (-dimensions => -1)', '-dimensions -1');
ok ($pso->setParams (-dimensions => 1), '-dimensions good');

for my $param (qw/numParticles/) {
    mustDie ("$pso->setParams (-$param => 0); $pso->init ()", "-$param zero");
    mustDie ("$pso->setParams (-$param => -1); $pso->init ()", "-$param neg");
    ok (($pso->setParams (-$param => 1), $pso->init ()), "-$param good");
}

for my $param (qw/inertia iterations meWeight numNeighbors stallSpeed themWeight/) {
    mustDie ("$pso->setParams (-$param => -1); $pso->init ()", "-$param neg");
    ok (($pso->setParams (-$param => 1), $pso->init ()), "-$param good");
}

mustDie ('$pso->setParams (-posMax => 0); $pso->setParams (-posMin => 0); $pso->init ()', '-posMax == -posMin');
mustDie ('$pso->setParams (-posMax => -1); $pso->setParams (-posMin => 0); $pso->init ()', '-posMax < -posMin');
ok ('$pso->setParams (-posMax => -1); $pso->setParams (-posMin => -2); $pso->init ()', '-posMax > -posMin');


sub fitFunc {
}


sub mustDie {
    my ( $test, $name ) = @_;

    eval $test;
    ok( defined $@, $name );
}
