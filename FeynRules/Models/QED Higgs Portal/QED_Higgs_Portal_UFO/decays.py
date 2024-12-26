# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.0.0 for Microsoft Windows (64-bit) (December 13, 2023)
# Date: Thu 26 Dec 2024 16:01:58


from object_library import all_decays, Decay
import particles as P


Decay_h = Decay(name = 'Decay_h',
                particle = P.h,
                partial_widths = {(P.chi__tilde__,P.chi):'(lambdachi**2*cmath.sqrt(-4*Mchi**2*Mh**2 + Mh**4))/(16.*cmath.pi*abs(Mh)**3)',
                                  (P.ele,P.ele__tilde__):'((-8*lambdah**2*ME**4 + 2*lambdah**2*ME**2*Mh**2)*cmath.sqrt(-4*ME**2*Mh**2 + Mh**4))/(16.*cmath.pi*abs(Mh)**3)',
                                  (P.mu,P.mu__tilde__):'((2*lambdah**2*Mh**2*MM**2 - 8*lambdah**2*MM**4)*cmath.sqrt(Mh**4 - 4*Mh**2*MM**2))/(16.*cmath.pi*abs(Mh)**3)',
                                  (P.ta,P.ta__tilde__):'((2*lambdah**2*Mh**2*ML**2 - 8*lambdah**2*ML**4)*cmath.sqrt(Mh**4 - 4*Mh**2*ML**2))/(16.*cmath.pi*abs(Mh)**3)'})

