import os
import sys
import subprocess
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm

one = map(float,"-0.000248206140000000	0.000508200490000000	0.00138963550000000	0.00178072280000000	0.00164791450000000	0.00164142170000000	0.00160329110000000	0.00108419860000000	0.000916052140000000	0.00124424230000000	0.00142714490000000	0.00195408500000000	0.00123125480000000	-0.000383033480000000	-0.00268420900000000	-0.00367596510000000	-0.00186788080000000	0.000671423510000000	0.00128258650000000	-0.000347320340000000	-0.00361541710000000	-0.00605389010000000	-0.00477329460000000	-0.00160721740000000	0.000697413870000000	0.00140198360000000	0.000509913780000000	-0.000772994940000000	-0.00141649970000000	-0.000696516130000000	0.000871995460000000	0.00146974940000000	-0.000127769540000000	0.000687437890000000	0.00115219220000000	0.00111815740000000	0.00114565190000000	0.00111193200000000	0.00115940240000000	0.00118258810000000	0.00154940910000000	-0.000278461490000000	-0.000322703330000000	0.000254559970000000	0.000564417450000000	-0.000324888590000000	-0.000459218370000000	-0.000380169720000000	-1.15323060000000e-05	0.000126507410000000	0.000439588880000000	-0.000436003610000000	-0.000476851190000000	-0.000478840140000000	-0.000385933530000000	-0.000180127990000000	0.000108142080000000	4.22443230000000e-05	-3.62348330000000e-05	3.28665560000000e-05	0.000105251330000000	5.69618280000000e-05	-0.000165074290000000	-0.000331975050000000	-0.000396901970000000	-0.000415214220000000	4.50077530000000e-05	0.000138373670000000	8.32978180000000e-05	-2.37624860000000e-05	-3.35429210000000e-05	-2.43879240000000e-05	-7.13937770000000e-05	5.11684410000000e-05	-9.17004650000000e-05	0.000112366670000000	0.000113622710000000	2.72696000000000e-05	-5.36118680000000e-05	1.80811580000000e-05	-3.44060430000000e-05	-0.000122471090000000	1.47911120000000e-06	0.000132086950000000	0.000127874270000000	-2.31706090000000e-05	-0.000264596430000000	-0.000107571770000000	0.000100138580000000	0.000132941770000000	2.53817600000000e-05	-3.40274780000000e-05	2.97698860000000e-05	-4.93223490000000e-05	-5.98815900000000e-05	7.03364930000000e-05	-0.000438228800000000	-0.000504878170000000	1.45866020000000e-05	0.000163365300000000	0.000103147690000000	-3.05184080000000e-05	-1.96171590000000e-05	-2.00941190000000e-05	8.28735720000000e-05	8.41752330000000e-05	-4.57439820000000e-06	-0.000315148240000000	-0.000387547830000000	-0.000386254160000000	-8.89222870000000e-05	0.000883655680000000	0.000731203710000000	-0.000492944090000000	-0.000522006940000000	-0.000530017890000000	-0.000294085300000000	6.96389400000000e-05	2.15545170000000e-05	-0.000387884880000000	-0.000570390260000000	-0.000417102740000000	0.000252876050000000	0.000369467540000000	0.00108105510000000	0.000679423680000000	-0.000195005340000000	3.26796290000000e-05".split('\t'))

two = [x for x in map(float,"6.0728e-02  -2.8465e-02   4.0315e-05  -7.1084e-05  -2.8210e-04  -1.5892e-04  -2.9249e-06   8.8681e-05   4.1803e-05  -3.7678e-04  -1.9632e-04   8.1989e-05  -5.8966e-05   8.9873e-05   5.1201e-05  -8.1240e-05   2.5059e-04  -1.6375e-04   2.6015e-04  -3.4187e-06   3.2514e-04  -4.9671e-04   3.4917e-05   3.0888e-05  -1.4810e-04  -9.8619e-05   2.1206e-04  -6.5177e-05  -1.3605e-05  -1.6820e-05   9.8165e-05   2.7446e-04   6.8879e-05  -5.3268e-05  -5.9838e-05   5.6147e-04   1.6358e-05   9.7536e-05  -1.6231e-04  -1.0849e-04  -3.0322e-04   1.1806e-04   1.0586e-04  -5.4383e-05  -1.6390e-08   3.4921e-06  -4.0464e-05  -8.2214e-05   2.0307e-04   3.0328e-05   1.7860e-04  -2.0808e-04   5.7529e-05   2.5681e-04   3.2688e-04  -1.7363e-04  -2.9062e-04  -7.3087e-05  -2.3758e-05  -1.9907e-06  -6.7034e-06   2.9695e-06   2.5441e-05   3.2610e-04  -6.1160e-04   2.0758e-04   1.3814e-04  -1.5335e-04   2.4614e-04  -1.7531e-04  -1.5214e-04   1.6134e-04   9.8156e-07   1.9119e-04  -3.1858e-04   1.7358e-04  -1.5364e-04  -1.4132e-04  -3.2287e-05   1.0957e-04   1.6327e-04  -4.8926e-06  -2.7889e-04  -1.7060e-05  -1.4486e-04   7.3797e-05  -4.7114e-05  -9.3366e-05   1.2303e-05  -3.4483e-05   3.1877e-04   1.6231e-04  -1.5561e-04  -1.3607e-04   3.4087e-04  -1.7361e-06  -2.6960e-04   3.6580e-04  -7.5748e-05  -5.0210e-05  -5.2966e-05  -9.4515e-05  -2.0064e-05   4.3477e-05  -4.5828e-05   1.5468e-05  -1.0921e-05   1.6420e-05  -1.2639e-05   8.2666e-05  -1.3898e-04  -1.2398e-04   1.6673e-05   1.2582e-04   1.3748e-04   7.4700e-05   8.4475e-06  -8.0269e-05   3.4473e-05   2.0991e-05  -5.2403e-05  -1.1068e-04  -6.6566e-05   5.5422e-05   6.0356e-05  -1.1841e-04   7.5837e-05   2.3525e-04".split('  '))]

three = [x for x in map(float,"3.9665e-04   1.9594e-02   2.8438e-07   1.0959e-04   1.6907e-04  -1.0609e-04  -5.9709e-05   4.8360e-05   1.0757e-04   1.0926e-04   4.7138e-05   4.5574e-06   6.9198e-05  -7.5191e-05   1.9739e-05   3.0644e-05  -6.1143e-06  -8.3707e-05  -6.7866e-05  -1.1197e-05  -1.5433e-04   9.4536e-05   2.5265e-04  -1.2052e-04  -1.0369e-04  -4.8149e-05   1.3169e-05  -3.5908e-06   5.2363e-05   1.6350e-05  -2.1275e-05  -2.7423e-05  -5.8823e-05   9.3754e-05  -1.2469e-04   1.4897e-04  -7.0698e-05   4.5566e-06   6.5150e-05   1.1797e-04   2.3494e-05   2.9939e-05  -8.7015e-05   3.4562e-05   5.9574e-05   6.4321e-05  -3.3367e-04   2.8843e-06   6.8625e-05  -2.3453e-05  -9.6210e-05   5.7793e-05   1.3354e-04  -2.5555e-04  -2.5228e-05  -4.0497e-05   9.3599e-05   5.3747e-05  -2.0394e-06  -3.4180e-05   5.6291e-06   1.7599e-06  -1.9361e-05  -1.3168e-04   2.8172e-04  -7.2662e-05  -2.1546e-04   1.5572e-04  -8.9376e-05   3.6020e-05   1.7179e-04  -1.6185e-04   4.9190e-06   3.9387e-05  -6.3531e-05  -1.0333e-05   2.0415e-05   5.8314e-05   6.1937e-05  -6.9411e-05   2.6733e-05   4.8274e-05   4.4397e-05  -6.9075e-05   3.6542e-05  -1.1973e-04   5.2467e-06   2.1024e-05  -3.4800e-05  -2.6039e-05   2.6672e-05  -1.1827e-06   5.7802e-05  -5.6001e-05  -4.7638e-05   2.4699e-05   1.8350e-05  -1.4780e-05   1.7256e-04  -5.2308e-05   1.0873e-04  -1.3695e-06  -1.0502e-05  -8.7420e-05   3.9930e-05   1.5259e-05   7.6885e-06   1.4524e-05  -1.6606e-06  -4.9175e-05   4.4502e-05  -9.9475e-05  -1.8451e-05  -1.1538e-04  -9.3272e-05   5.4999e-05   2.6652e-05   4.5891e-06  -6.7009e-06  -1.5656e-06  -7.3672e-06  -1.1709e-05   6.6960e-05   2.4801e-05   4.2509e-05  -8.9745e-05   8.1479e-06  -1.6059e-04".split('  '))]

n_two = two[96:] + two[:96]
n_three = three[96:] + three[:96]

plt.plot([x for x in range(len(one))],one, hold=True)
plt.plot([x for x in range(len(one))], two, hold=True)
plt.plot([x for x in range(len(one))], three)
plt.show()
