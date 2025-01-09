__date__ = "30 December 2024"
__author__ = "hypercube256@snu.ac.kr"

import cmath
from object_library import all_functions, Function

#
# shortcuts for functions from cmath
#

complexconjugate = Function(name = 'complexconjugate',
                            arguments = ('z',),
                            expression = 'z.conjugate()')


re = Function(name = 're',
              arguments = ('z',),
              expression = 'z.real')

im = Function(name = 'im',
              arguments = ('z',),
              expression = 'z.imag')

# New functions (trigonometric)

sec = Function(name = 'sec',
             arguments = ('z',),
             expression = '1./cmath.cos(z.real)')

asec = Function(name = 'asec',
             arguments = ('z',),
             expression = 'cmath.acos(1./(z.real))')

csc = Function(name = 'csc',
             arguments = ('z',),
             expression = '1./cmath.sin(z.real)')

acsc = Function(name = 'acsc',
             arguments = ('z',),
             expression = 'cmath.asin(1./(z.real))')

cot = Function(name = 'cot',
               arguments = ('z',),
               expression = '1./cmath.tan(z.real)')

# Heaviside theta function

theta_function = Function(name = 'theta_function',
             arguments = ('x','y','z'),
             expression = 'y if x else z')

# Auxiliary functions for NLO

scalarF = Function(name = 'scalarF',
            arguments = ('m','s'),
            expression = '2. - 2.*cmath.log(m) - (2.*cmath.sqrt(4.*m**2-s)/cmath.sqrt(s))*cmath.atan(cmath.sqrt(s)/cmath.sqrt(4.*m**2-s))'
)

scalarFP = Function(name = 'scalarFP',
            arguments = ('m','s'),
            expression = '-(1./s) + (4.*m**2/(cmath.sqrt(4.*m**2-s)*cmath.sqrt(s**3)))*(cmath.atan(cmath.sqrt(s)/cmath.sqrt(4.*m**2-s)))'
)
