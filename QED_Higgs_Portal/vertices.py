# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.0.0 for Microsoft Windows (64-bit) (December 13, 2023)
# Date: Mon 30 Dec 2024 02:28:48


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.e__tilde__, P.e, P.A ],
             color = [ '1' ],
             lorentz = [ L.FFV1 ],
             couplings = {(0,0):C.GC_1})

V_2 = Vertex(name = 'V_2',
             particles = [ P.m__tilde__, P.m, P.A ],
             color = [ '1' ],
             lorentz = [ L.FFV1 ],
             couplings = {(0,0):C.GC_1})

V_3 = Vertex(name = 'V_3',
             particles = [ P.tt__tilde__, P.tt, P.A ],
             color = [ '1' ],
             lorentz = [ L.FFV1 ],
             couplings = {(0,0):C.GC_1})

V_4 = Vertex(name = 'V_4',
             particles = [ P.e__tilde__, P.e, P.H ],
             color = [ '1' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_2})

V_5 = Vertex(name = 'V_5',
             particles = [ P.m__tilde__, P.m, P.H ],
             color = [ '1' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_3})

V_6 = Vertex(name = 'V_6',
             particles = [ P.tt__tilde__, P.tt, P.H ],
             color = [ '1' ],
             lorentz = [ L.FFS1 ],
             couplings = {(0,0):C.GC_4})

