import matrix2SomeNamepy as matrix2py
def invert_momenta(p):
   """ fortran/C-python do not order table in the same order"""
   new_p = []
   for i in range(len(p[0])):  new_p.append([0]*len(p))
   for i, onep in enumerate(p):
      for j, x in enumerate(onep):
         new_p[j][i] = x
   return new_p
matrix2py.initialise('param_card.dat')
nhel = 0 # sum over all helicities
# 2HDM setup: mX=600, sba=1, tanb=0.3, type=2 --> topo=P0_gg_ttx, mtt=506.969, me2XX/me2SM=0.516962/0.517445=0.999065
# ME ids are:  [21, 21, 6, -6]
# Q = 245.40846875
alphaS = 0.112812567332
p=[
  [501.564,0,0,501.564],
  [128.108,0,0,-128.108],
  [360.934,50.0707,167.275,265.174],
  [268.738,-50.0707,-167.275,108.282],
]
P = invert_momenta(p)
me2 = matrix2py.get_me(P,alphaS,nhel)
print "ME2(ttx)=",me2

# ### MEcode=0, topology=ttx, Q=207.957183669, alphaS=0.10498082074 ->
# p = [[148.7495, 0.0, 0.0, 148.7495], [346.258375, 0.0, 0.0, -346.258375],
#      [285.4675, -99.5479765625, 60.29248046875, -197.2335625], [209.539609375, 99.5479765625, -60.29248046875, -0.2753184814453125]]
# alphas=0.10498082074
# P = invert_momenta(p)
# print "p =",p
# print "P =",P
# me2 = matrix2py.get_me(P,alphas,nhel)
# print "ME2(ttx)=",me2


### MEcode=1, topology=ttxg, Q=209.048729925, alphaS=0.104907446345 ->
# p = [[264.3046875, 0.0, 0.0, 264.3046875], [328.02553125, 0.0, 0.0, -328.02553125],
#     [242.8330625, 18.229767578125, 58.28096484375, 157.198890625], [199.76946875, -59.23303515625, -29.940345703125, -79.73228125], [149.727625, 41.00326953125, -28.34062109375, -141.18746875]]
# alphas = 0.104907446345                                                                                                       
# P = invert_momenta(p)
# me2 = matrix2py.get_me(P,alphas,nhel)
# print "ME2(ttxg)=",me2


# ### MEcode=2, topology=ttxgg, Q=268.652716997, alphaS=0.101510739152 ->
# p = [[596.8451875, 0.0, 0.0, 596.8451875], [245.78221875, 0.0, 0.0, -245.78221875],
#      [433.40378125, -140.238125, -88.751625, 361.48684375], [238.4015625, 154.838421875, 24.286923828125, 54.7362890625],
#      [56.67294140625, -6.11514990234375, 34.39092578125, 44.62837109375], [114.1488359375, -8.485150390625, 30.0737734375, -109.7885546875]]
# alphas = 0.101510739152
# P = invert_momenta(p)
# print "P =",P
# me2 = matrix2py.get_me(P,alphas,nhel)
# print "ME2(ttxgg)=",me2
