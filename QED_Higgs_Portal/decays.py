# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.0.0 for Microsoft Windows (64-bit) (December 13, 2023)
# Date: Mon 30 Dec 2024 02:28:48


from object_library import all_decays, Decay
import particles as P


Decay_H = Decay(name = 'Decay_H',
                particle = P.H,
                partial_widths = {(P.e,P.e__tilde__):'((-8*lambdah**2*Me**4 + 2*lambdah**2*Me**2*MH**2)*cmath.sqrt(-4*Me**2*MH**2 + MH**4))/(16.*cmath.pi*abs(MH)**3)',
                                  (P.m,P.m__tilde__):'((2*lambdah**2*MH**2*MM**2 - 8*lambdah**2*MM**4)*cmath.sqrt(MH**4 - 4*MH**2*MM**2))/(16.*cmath.pi*abs(MH)**3)',
                                  (P.tt,P.tt__tilde__):'((2*lambdah**2*MH**2*MTA**2 - 8*lambdah**2*MTA**4)*cmath.sqrt(MH**4 - 4*MH**2*MTA**2))/(16.*cmath.pi*abs(MH)**3)'})

