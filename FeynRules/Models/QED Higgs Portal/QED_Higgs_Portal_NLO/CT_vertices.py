# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.0.0 for Microsoft Windows (64-bit) (December 13, 2023)
# Date: Sat 28 Dec 2024 18:41:08


from object_library import all_vertices, all_CTvertices, Vertex, CTVertex
import particles as P
import CT_couplings as C
import lorentz as L


V_1 = CTVertex(name = 'V_1',
               type = 'R2',
               particles = [ P.ele__tilde__, P.ele, P.A ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               loop_particles = [ [ [P.A, P.ele] ], [ [P.ele, P.h] ] ],
               couplings = {(0,0,0):C.R2GC_17_7,(0,0,1):C.R2GC_17_8})

V_2 = CTVertex(name = 'V_2',
               type = 'R2',
               particles = [ P.mu__tilde__, P.mu, P.A ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               loop_particles = [ [ [P.A, P.mu] ], [ [P.h, P.mu] ] ],
               couplings = {(0,0,0):C.R2GC_17_7,(0,0,1):C.R2GC_18_9})

V_3 = CTVertex(name = 'V_3',
               type = 'R2',
               particles = [ P.ta__tilde__, P.ta, P.A ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               loop_particles = [ [ [P.A, P.ta] ], [ [P.h, P.ta] ] ],
               couplings = {(0,0,0):C.R2GC_17_7,(0,0,1):C.R2GC_19_10})

V_4 = CTVertex(name = 'V_4',
               type = 'R2',
               particles = [ P.ele__tilde__, P.ele, P.h ],
               color = [ '1' ],
               lorentz = [ L.FFS1 ],
               loop_particles = [ [ [P.A, P.ele] ], [ [P.ele, P.h] ] ],
               couplings = {(0,0,0):C.R2GC_32_18,(0,0,1):C.R2GC_32_19})

V_5 = CTVertex(name = 'V_5',
               type = 'R2',
               particles = [ P.mu__tilde__, P.mu, P.h ],
               color = [ '1' ],
               lorentz = [ L.FFS1 ],
               loop_particles = [ [ [P.A, P.mu] ], [ [P.h, P.mu] ] ],
               couplings = {(0,0,0):C.R2GC_33_20,(0,0,1):C.R2GC_33_21})

V_6 = CTVertex(name = 'V_6',
               type = 'R2',
               particles = [ P.ta__tilde__, P.ta, P.h ],
               color = [ '1' ],
               lorentz = [ L.FFS1 ],
               loop_particles = [ [ [P.A, P.ta] ], [ [P.h, P.ta] ] ],
               couplings = {(0,0,0):C.R2GC_38_26,(0,0,1):C.R2GC_38_27})

V_7 = CTVertex(name = 'V_7',
               type = 'R2',
               particles = [ P.ele__tilde__, P.ele ],
               color = [ '1' ],
               lorentz = [ L.FF2, L.FF3 ],
               loop_particles = [ [ [P.A, P.ele] ] ],
               couplings = {(0,0,0):C.R2GC_23_12,(0,1,0):C.R2GC_22_11})

V_8 = CTVertex(name = 'V_8',
               type = 'R2',
               particles = [ P.mu__tilde__, P.mu ],
               color = [ '1' ],
               lorentz = [ L.FF2, L.FF3 ],
               loop_particles = [ [ [P.A, P.mu] ] ],
               couplings = {(0,0,0):C.R2GC_25_13,(0,1,0):C.R2GC_22_11})

V_9 = CTVertex(name = 'V_9',
               type = 'R2',
               particles = [ P.ta__tilde__, P.ta ],
               color = [ '1' ],
               lorentz = [ L.FF2, L.FF3 ],
               loop_particles = [ [ [P.A, P.ta] ] ],
               couplings = {(0,0,0):C.R2GC_37_25,(0,1,0):C.R2GC_22_11})

V_10 = CTVertex(name = 'V_10',
                type = 'R2',
                particles = [ P.h, P.h ],
                color = [ '1' ],
                lorentz = [ L.SS1, L.SS3 ],
                loop_particles = [ [ [P.ele] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,0):C.R2GC_34_22,(0,0,1):C.R2GC_34_23,(0,0,2):C.R2GC_34_24,(0,1,0):C.R2GC_30_15,(0,1,1):C.R2GC_30_16,(0,1,2):C.R2GC_30_17})

V_11 = CTVertex(name = 'V_11',
                type = 'R2',
                particles = [ P.A, P.A ],
                color = [ '1' ],
                lorentz = [ L.VV2, L.VV3 ],
                loop_particles = [ [ [P.ele] ], [ [P.ele], [P.mu], [P.ta] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,0):C.R2GC_7_29,(0,0,2):C.R2GC_7_30,(0,0,3):C.R2GC_7_31,(0,1,1):C.R2GC_26_14})

V_12 = CTVertex(name = 'V_12',
                type = 'R2',
                particles = [ P.A, P.A, P.h ],
                color = [ '1' ],
                lorentz = [ L.VVS1 ],
                loop_particles = [ [ [P.ele] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,0):C.R2GC_8_32,(0,0,1):C.R2GC_8_33,(0,0,2):C.R2GC_8_34})

V_13 = CTVertex(name = 'V_13',
                type = 'R2',
                particles = [ P.h, P.h, P.h ],
                color = [ '1' ],
                lorentz = [ L.SSS1 ],
                loop_particles = [ [ [P.ele] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,0):C.R2GC_15_1,(0,0,1):C.R2GC_15_2,(0,0,2):C.R2GC_15_3})

V_14 = CTVertex(name = 'V_14',
                type = 'R2',
                particles = [ P.A, P.A, P.A, P.A ],
                color = [ '1' ],
                lorentz = [ L.VVVV1 ],
                loop_particles = [ [ [P.ele], [P.mu], [P.ta] ] ],
                couplings = {(0,0,0):C.R2GC_6_28})

V_15 = CTVertex(name = 'V_15',
                type = 'R2',
                particles = [ P.A, P.A, P.h, P.h ],
                color = [ '1' ],
                lorentz = [ L.VVSS1 ],
                loop_particles = [ [ [P.ele] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,0):C.R2GC_9_35,(0,0,1):C.R2GC_9_36,(0,0,2):C.R2GC_9_37})

V_16 = CTVertex(name = 'V_16',
                type = 'R2',
                particles = [ P.h, P.h, P.h, P.h ],
                color = [ '1' ],
                lorentz = [ L.SSSS1 ],
                loop_particles = [ [ [P.ele] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,0):C.R2GC_16_4,(0,0,1):C.R2GC_16_5,(0,0,2):C.R2GC_16_6})

V_17 = CTVertex(name = 'V_17',
                type = 'UV',
                particles = [ P.cc__tilde__, P.cc, P.h ],
                color = [ '1' ],
                lorentz = [ L.SSS1 ],
                loop_particles = [ [ [P.cc] ], [ [P.cc, P.h] ], [ [P.ele] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,0):C.UVGC_31_41,(0,0,2):C.UVGC_31_42,(0,0,3):C.UVGC_31_43,(0,0,4):C.UVGC_31_44,(0,0,1):C.UVGC_31_45})

V_18 = CTVertex(name = 'V_18',
                type = 'UV',
                particles = [ P.ele__tilde__, P.ele, P.A ],
                color = [ '1' ],
                lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.A, P.ele] ], [ [P.ele] ], [ [P.ele, P.h] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,0):C.UVGC_11_2,(0,0,2):C.UVGC_17_11,(0,1,1):C.UVGC_28_30,(0,1,3):C.UVGC_28_31,(0,1,4):C.UVGC_28_32,(0,1,0):C.UVGC_28_33,(0,1,2):C.UVGC_28_34,(0,2,1):C.UVGC_28_30,(0,2,3):C.UVGC_28_31,(0,2,4):C.UVGC_28_32,(0,2,0):C.UVGC_28_33,(0,2,2):C.UVGC_28_34})

V_19 = CTVertex(name = 'V_19',
                type = 'UV',
                particles = [ P.mu__tilde__, P.mu, P.A ],
                color = [ '1' ],
                lorentz = [ L.FFV1, L.FFV4, L.FFV5 ],
                loop_particles = [ [ [P.A, P.mu] ], [ [P.ele] ], [ [P.h, P.mu] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,2):C.UVGC_18_12,(0,2,0):C.UVGC_11_2,(0,1,1):C.UVGC_28_30,(0,1,3):C.UVGC_28_31,(0,1,4):C.UVGC_28_32,(0,1,0):C.UVGC_29_35,(0,1,2):C.UVGC_29_36})

V_20 = CTVertex(name = 'V_20',
                type = 'UV',
                particles = [ P.ta__tilde__, P.ta, P.A ],
                color = [ '1' ],
                lorentz = [ L.FFV1, L.FFV4, L.FFV5 ],
                loop_particles = [ [ [P.A, P.ta] ], [ [P.ele] ], [ [P.h, P.ta] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,2):C.UVGC_19_13,(0,2,0):C.UVGC_11_2,(0,1,1):C.UVGC_28_30,(0,1,3):C.UVGC_28_31,(0,1,4):C.UVGC_28_32,(0,1,0):C.UVGC_36_64,(0,1,2):C.UVGC_36_65})

V_21 = CTVertex(name = 'V_21',
                type = 'UV',
                particles = [ P.ele__tilde__, P.ele, P.h ],
                color = [ '1' ],
                lorentz = [ L.FFS1 ],
                loop_particles = [ [ [P.A, P.ele] ], [ [P.cc] ], [ [P.ele] ], [ [P.ele, P.h] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,1):C.UVGC_32_46,(0,0,2):C.UVGC_32_47,(0,0,4):C.UVGC_32_48,(0,0,5):C.UVGC_32_49,(0,0,0):C.UVGC_32_50,(0,0,3):C.UVGC_32_51})

V_22 = CTVertex(name = 'V_22',
                type = 'UV',
                particles = [ P.mu__tilde__, P.mu, P.h ],
                color = [ '1' ],
                lorentz = [ L.FFS1 ],
                loop_particles = [ [ [P.A, P.mu] ], [ [P.cc] ], [ [P.ele] ], [ [P.h, P.mu] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,1):C.UVGC_33_52,(0,0,2):C.UVGC_33_53,(0,0,4):C.UVGC_33_54,(0,0,5):C.UVGC_33_55,(0,0,0):C.UVGC_33_56,(0,0,3):C.UVGC_33_57})

V_23 = CTVertex(name = 'V_23',
                type = 'UV',
                particles = [ P.ta__tilde__, P.ta, P.h ],
                color = [ '1' ],
                lorentz = [ L.FFS1 ],
                loop_particles = [ [ [P.A, P.ta] ], [ [P.cc] ], [ [P.ele] ], [ [P.h, P.ta] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,1):C.UVGC_38_68,(0,0,2):C.UVGC_38_69,(0,0,4):C.UVGC_38_70,(0,0,5):C.UVGC_38_71,(0,0,0):C.UVGC_38_72,(0,0,3):C.UVGC_38_73})

V_24 = CTVertex(name = 'V_24',
                type = 'UV',
                particles = [ P.ele__tilde__, P.ele ],
                color = [ '1' ],
                lorentz = [ L.FF1, L.FF2, L.FF3 ],
                loop_particles = [ [ [P.A, P.ele] ], [ [P.ele, P.h] ] ],
                couplings = {(0,1,0):C.UVGC_23_18,(0,1,1):C.UVGC_23_19,(0,2,0):C.UVGC_22_16,(0,2,1):C.UVGC_22_17,(0,0,1):C.UVGC_10_1})

V_25 = CTVertex(name = 'V_25',
                type = 'UV',
                particles = [ P.mu__tilde__, P.mu ],
                color = [ '1' ],
                lorentz = [ L.FF1, L.FF2, L.FF3 ],
                loop_particles = [ [ [P.A, P.mu] ], [ [P.h, P.mu] ] ],
                couplings = {(0,1,0):C.UVGC_25_22,(0,1,1):C.UVGC_25_23,(0,2,0):C.UVGC_24_20,(0,2,1):C.UVGC_24_21,(0,0,1):C.UVGC_12_3})

V_26 = CTVertex(name = 'V_26',
                type = 'UV',
                particles = [ P.ta__tilde__, P.ta ],
                color = [ '1' ],
                lorentz = [ L.FF1, L.FF2, L.FF3 ],
                loop_particles = [ [ [P.A, P.ta] ], [ [P.h, P.ta] ] ],
                couplings = {(0,1,0):C.UVGC_37_66,(0,1,1):C.UVGC_37_67,(0,2,0):C.UVGC_35_62,(0,2,1):C.UVGC_35_63,(0,0,1):C.UVGC_14_4})

V_27 = CTVertex(name = 'V_27',
                type = 'UV',
                particles = [ P.h, P.h ],
                color = [ '1' ],
                lorentz = [ L.SS1, L.SS3 ],
                loop_particles = [ [ [P.cc] ], [ [P.ele] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,0):C.UVGC_34_58,(0,0,1):C.UVGC_34_59,(0,0,2):C.UVGC_34_60,(0,0,3):C.UVGC_34_61,(0,1,0):C.UVGC_30_37,(0,1,1):C.UVGC_30_38,(0,1,2):C.UVGC_30_39,(0,1,3):C.UVGC_30_40})

V_28 = CTVertex(name = 'V_28',
                type = 'UV',
                particles = [ P.A, P.A ],
                color = [ '1' ],
                lorentz = [ L.VV1, L.VV3 ],
                loop_particles = [ [ [P.ele] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,1,0):C.UVGC_26_24,(0,1,1):C.UVGC_26_25,(0,1,2):C.UVGC_26_26,(0,0,0):C.UVGC_27_27,(0,0,1):C.UVGC_27_28,(0,0,2):C.UVGC_27_29})

V_29 = CTVertex(name = 'V_29',
                type = 'UV',
                particles = [ P.h, P.h, P.h ],
                color = [ '1' ],
                lorentz = [ L.SSS1 ],
                loop_particles = [ [ [P.ele] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,0):C.UVGC_15_5,(0,0,1):C.UVGC_15_6,(0,0,2):C.UVGC_15_7})

V_30 = CTVertex(name = 'V_30',
                type = 'UV',
                particles = [ P.h, P.h, P.h, P.h ],
                color = [ '1' ],
                lorentz = [ L.SSSS1 ],
                loop_particles = [ [ [P.ele] ], [ [P.mu] ], [ [P.ta] ] ],
                couplings = {(0,0,0):C.UVGC_16_8,(0,0,1):C.UVGC_16_9,(0,0,2):C.UVGC_16_10})

V_31 = CTVertex(name = 'V_31',
                type = 'UV',
                particles = [ P.cc__tilde__, P.cc ],
                color = [ '1' ],
                lorentz = [ L.SS1, L.SS2 ],
                loop_particles = [ [ [P.cc, P.h] ] ],
                couplings = {(0,0,0):C.UVGC_21_15,(0,1,0):C.UVGC_20_14})

