from object_library import all_decays, Decay
import particles as P


Decay_H = Decay(name = 'Decay_H',
                particle = P.H,
                partial_widths = {(P.e__minus__,P.e__plus__):'((-8*lambdah**2*Me**4 + 2*lambdah**2*Me**2*MH**2)*cmath.sqrt(-4*Me**2*MH**2 + MH**4))/(16.*cmath.pi*abs(MH)**3)',
                                  (P.m__minus__,P.m__plus__):'((2*lambdah**2*MH**2*MM**2 - 8*lambdah**2*MM**4)*cmath.sqrt(MH**4 - 4*MH**2*MM**2))/(16.*cmath.pi*abs(MH)**3)',
                                  (P.tt__minus__,P.tt__plus__):'((2*lambdah**2*MH**2*MTA**2 - 8*lambdah**2*MTA**4)*cmath.sqrt(MH**4 - 4*MH**2*MTA**2))/(16.*cmath.pi*abs(MH)**3)'})

