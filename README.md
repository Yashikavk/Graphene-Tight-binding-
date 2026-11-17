# Graphene-Tight-binding
If an electron sits in Graphene, what energies can it have and how does it vary with momentum? The tight-binding model calculates this exact relationship, revealing a gapless band structure where electrons behave like massless particles near the Dirac points. 

Graphene is a honeycomb lattice with two carbon atoms (A and B) per unit cell.
One pz orbital per atom, electrons hop only between nearest neighbours with
energy -t. In the basis (|A>, |B>) the Bloch Hamiltonian is

        H(k) = [[ delta,      -t f(k) ],
                [ -t f*(k),   -delta  ]],      f(k) = sum_j exp(i k . d_j)

where d_j are the three vectors from an A atom to its B neighbours.
For delta = 0:  E(k) = +/- t |f(k)|.  
Setting delta = 0 breaks the A/B symmetry and opens a gap.
