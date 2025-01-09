# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.0.0 for Microsoft Windows (64-bit) (December 13, 2023)
# Date: Sat 28 Dec 2024 18:41:08


from object_library import all_decays, Decay
import particles as P


Decay_h = Decay(name = 'Decay_h',
                particle = P.h,
                partial_widths = {(P.cc__tilde__,P.cc):'(lambdacc**2*cmath.sqrt(-4*Mcc**2*Mh**2 + Mh**4))/(16.*cmath.pi*abs(Mh)**3)',
                                  (P.ele,P.ele__tilde__):'((-8*lambdah**2*Mele**4 + 2*lambdah**2*Mele**2*Mh**2)*cmath.sqrt(-4*Mele**2*Mh**2 + Mh**4))/(16.*cmath.pi*abs(Mh)**3)',
                                  (P.mu,P.mu__tilde__):'((2*lambdah**2*Mh**2*Mmu**2 - 8*lambdah**2*Mmu**4)*cmath.sqrt(Mh**4 - 4*Mh**2*Mmu**2))/(16.*cmath.pi*abs(Mh)**3)',
                                  (P.ta,P.ta__tilde__):'((2*lambdah**2*Mh**2*Mtau**2 - 8*lambdah**2*Mtau**4)*cmath.sqrt(Mh**4 - 4*Mh**2*Mtau**2))/(16.*cmath.pi*abs(Mh)**3)'})

