# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.0.0 for Microsoft Windows (64-bit) (December 13, 2023)
# Date: Fri 20 Dec 2024 11:45:40


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.f__tilde__, P.f, P.h ],
             color = [ '1' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_2})

V_2 = Vertex(name = 'V_2',
             particles = [ P.ele__tilde__, P.ele, P.h ],
             color = [ '1' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_3})

V_3 = Vertex(name = 'V_3',
             particles = [ P.mu__tilde__, P.mu, P.h ],
             color = [ '1' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_3})

V_4 = Vertex(name = 'V_4',
             particles = [ P.ta__tilde__, P.ta, P.h ],
             color = [ '1' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_3})

V_5 = Vertex(name = 'V_5',
             particles = [ P.ele__tilde__, P.ele, P.A ],
             color = [ '1' ],
             lorentz = [ L.FFV1 ],
             couplings = {(0,0):C.GC_1})

V_6 = Vertex(name = 'V_6',
             particles = [ P.mu__tilde__, P.mu, P.A ],
             color = [ '1' ],
             lorentz = [ L.FFV1 ],
             couplings = {(0,0):C.GC_1})

V_7 = Vertex(name = 'V_7',
             particles = [ P.ta__tilde__, P.ta, P.A ],
             color = [ '1' ],
             lorentz = [ L.FFV1 ],
             couplings = {(0,0):C.GC_1})

