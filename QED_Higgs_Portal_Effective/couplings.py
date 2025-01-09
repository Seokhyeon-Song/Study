from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot



GC_1 = Coupling(name = 'GC_1',
                value = 'ee*complex(0,1)',
                order = {'QED':1})

GC_2 = Coupling(name = 'GC_2',
                value = '-(complex(0,1)*lambdah*Me)',
                order = {'QED':1})

GC_3 = Coupling(name = 'GC_3',
                value = '-(complex(0,1)*lambdah*MM)',
                order = {'QED':1})

GC_4 = Coupling(name = 'GC_4',
                value = '-(complex(0,1)*lambdah*MTA)',
                order = {'QED':1})

