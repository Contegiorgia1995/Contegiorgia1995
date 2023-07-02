# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 17:38:08 2022

@author: Giorgia
"""

ife = 5
educ = 3
INNO_pos = range(0,5)

Spgrid = array([-4.49557054e-01, -2.29935939e-01, -9.68375446e-02, -2.85940025e-02,
       -3.53744357e-03,  1.16850017e-09,  3.68620065e-03,  2.91890238e-02,
        9.81763396e-02,  2.32316017e-01,  4.53275925e-01,  7.82723932e-01,
        1.24232791e+00,  1.85375572e+00,  2.63867524e+00,  3.61875434e+00,
        4.81566087e+00,  6.25106273e+00,  7.94662776e+00,  9.92402385e+00,
        1.22049189e+01,  1.48109806e+01,  1.77638771e+01,  2.10852761e+01,
        2.47968455e+01,  2.89202531e+01,  3.34771669e+01,  3.84892546e+01,
        4.39781843e+01,  4.99656237e+01,  5.64732407e+01,  6.35227031e+01,
        7.11356790e+01,  7.93338360e+01,  8.81388421e+01,  9.75723652e+01,
        1.07656073e+02,  1.18411634e+02,  1.29860715e+02,  1.42024984e+02,
        1.54926110e+02,  1.68585760e+02,  1.83025602e+02,  1.98267304e+02,
        2.14332534e+02,  2.31242960e+02,  2.49020249e+02,  2.67686069e+02,
        2.87262088e+02,  3.07769975e+02,  3.29231397e+02,  3.51668021e+02,
        3.75101517e+02,  3.99553550e+02,  4.25045791e+02,  4.51599905e+02,
        4.79237562e+02,  5.07980429e+02,  5.37850174e+02,  5.68868465e+02,
        6.01056969e+02,  6.34437355e+02,  6.69031290e+02,  7.04860443e+02,
        7.41946481e+02,  7.80311071e+02,  8.19975883e+02,  8.60962583e+02,
        9.03292840e+02,  9.46988322e+02,  9.92070696e+02,  1.03856163e+03,
        1.08648279e+03,  1.13585585e+03,  1.18670247e+03,  1.23904433e+03,
        1.29290308e+03,  1.34830040e+03,  1.40525796e+03,  1.46379742e+03])

values = array([[ 10.47912691,  11.0922498 ,  11.34450658,  11.45837614,
         11.49818626,  11.50372841,  11.50851153,  11.54120272,
         11.62822745,  11.79526906,  12.06505511,  12.45408715,
         12.97635511,  13.64262681,  14.45479056,  15.40588675,
         16.49238573,  17.71545221,  19.06716317,  20.53047931,
         22.08970077,  23.75607508,  25.52371297,  27.38319239,
         29.29596822,  31.28677169,  33.35740207,  35.50476163,
         37.69942109,  39.9561176 ,  42.28090377,  44.67235486,
         47.1144836 ,  49.60202688,  52.1423648 ,  54.73890734,
         57.38485797,  60.07053646,  62.79898439,  65.58395064,
         68.42299356,  71.3071313 ,  74.21849637,  38.92316483,
         40.0736902 ,  41.23919056,  42.41946973,  43.6143394 ,
         44.82361858,  46.04713309,  47.28471519,  48.53620312,
         49.80144074,  51.08027725,  52.37256683,  53.67816839,
         54.9969453 ,  56.32876516,  57.67349956,  59.03102392,
         60.40121726,  61.78396204,  63.17914398,  64.58665194,
         66.00637776,  67.43821611,  68.88206439,  70.3378226 ,
         71.80539323,  73.28468117,  74.77559358,  76.27803985,
         77.79193145,  79.31718189,  80.85370666,  82.40142308,
         83.96025032,  85.53010928,  87.11092254,  88.70261431],
       [ 11.46060136,  11.79436429,  11.99046942,  12.09007323,
         12.12661602,  12.13177387,  12.13616197,  12.16708576,
         12.24799222,  12.40479245,  12.66121475,  13.0308613 ,
         13.53236635,  14.1762262 ,  14.96457577,  15.8798035 ,
         16.9360332 ,  18.12992836,  19.45179962,  20.88405396,
         22.41909822,  24.06314767,  25.80998018,  27.64816829,
         29.54480628,  31.52277265,  33.57765564,  35.71057482,
         37.89246785,  40.13745091,  42.45132978,  44.83269386,
         47.26682589,  49.74684816,  52.27959685,  54.8691674 ,
         57.50987843,  60.19121667,  62.91361498,  65.69285366,
         68.52731081,  71.4066754 ,  74.30819907,  77.26016757,
         80.2614557 ,  83.31027443,  86.38300099,  89.49445841,
         92.6536986 ,  95.85970642,  99.09992665, 102.3693898 ,
        105.67954898, 109.03386139, 112.42715032, 115.8466127 ,
        119.29665266, 122.78992582, 126.32453819, 129.88933744,
        133.47751466, 137.10586996, 140.77469672, 144.48046262,
        148.20180189, 151.94989086, 155.73585001, 159.56058596,
        163.41449834, 167.2924661 , 171.20668713, 175.15902727,
        179.14572641, 183.15402079, 187.19011737, 191.2625006 ,
        195.36646217, 199.49613704, 203.64556936, 207.82706532],
       [ 12.57540977,  12.87430012,  13.05366673,  13.14549822,
         13.17919295,  13.18394896,  13.18806619,  13.21655414,
         13.29112843,  13.43522244,  13.67196871,  14.02284918,
         14.48665725,  15.09064985,  15.8237386 ,  16.69452456,
         17.70955431,  18.86448385,  20.14417943,  21.51888441,
         23.01084796,  24.61527687,  26.32500372,  28.12075782,
         29.9873513 ,  31.9416024 ,  33.979487  ,  36.09564303,
         38.25022506,  40.4744763 ,  42.76997976,  45.134377  ,
         47.55063845,  50.01517952,  52.53466107,  55.11190197,
         57.74033933,  60.41024919,  63.12328196,  65.89373422,
         68.71964275,  71.59345408,  74.48783501,  77.43262455,
         80.42711777,  83.46950623,  86.53554243,  89.64142304,
         92.79537594,  95.99635497,  99.23160239, 102.4965479 ,
        105.8024552 , 109.15271007, 112.54216087, 115.95801874,
        119.40457522, 122.89450447, 126.42506984, 129.98674217,
        133.57210565, 137.19750747, 140.86350712, 144.5666263 ,
        148.28605607, 152.03170732, 155.81528251, 159.63767877,
        163.48994793, 167.36584417, 171.27796027, 175.2282661 ,
        179.21355283, 183.22029207, 187.25456387, 191.32256353,
        195.42470344, 199.55304133, 203.70146924, 207.8813768 ],
       [ 14.36611927,  14.63212092,  14.79299227,  14.87538102,
         14.90561535,  14.9098831 ,  14.9135757 ,  14.93911785,
         15.00400551,  15.12988516,  15.33716147,  15.64477477,
         16.07096665,  16.60605889,  17.287982  ,  18.08895199,
         19.03696757,  20.1157281 ,  21.2992791 ,  22.61044515,
         24.04569586,  25.59730671,  27.23998359,  28.96188332,
         30.77404733,  32.67856768,  34.67097392,  36.74217473,
         38.86785269,  41.05947471,  43.32051423,  45.65298845,
         48.038552  ,  50.47570636,  52.9703128 ,  55.52444616,
         58.13173097,  60.78245546,  63.47684255,  66.22967334,
         69.03898536,  71.88965892,  74.76958749,  77.69954136,
         80.68006674,  83.70932249,  86.76782259,  89.86164865,
         93.00400792,  96.19396185,  99.42342608, 102.67975087,
        105.97597428, 109.3170553 , 112.70106346, 116.10460371,
        119.53798688, 123.01377908, 126.53286055, 130.08855656,
        133.67143511, 137.28885208, 140.94737917, 144.6444155 ,
        148.36737241, 152.10681867, 155.88384967, 159.69935034,
        163.55273326, 167.42237206, 171.3228613 , 175.26170833,
        179.23831688, 183.23885824, 187.26535196, 191.3253659 ,
        195.42150307, 199.54824602, 203.69796184, 207.87270021],
       [ 17.14290005,  17.36593977,  17.49736861,  17.56471251,
         17.58943142,  17.59292088,  17.5959436 ,  17.61685383,
         17.67339759,  17.78325843,  17.96398984,  18.23293101,
         18.60591036,  19.07407826,  19.66732199,  20.39037054,
         21.21865875,  22.18597165,  23.27458582,  24.48450126,
         25.82642569,  27.2762954 ,  28.80682103,  30.44265579,
         32.18302305,  34.02213668,  35.9349472 ,  37.91605043,
         39.97453956,  42.11357927,  44.33128734,  46.6212353 ,
         48.96462171,  51.35201166,  53.80402935,  56.31853056,
         58.88534766,  61.49808504,  64.1628316 ,  66.88788759,
         69.670097  ,  72.5006388 ,  75.3629608 ,  78.27388084,
         81.23297696,  84.24188302,  87.27600229,  90.35247294,
         93.47844012,  96.65277957,  99.86411532, 103.10569661,
        106.3888394 , 109.71738123, 113.08815106, 116.48609301,
        119.90973886, 123.37510291, 126.88424063, 130.42959382,
        134.00161479, 137.61045535, 141.26066358, 144.94940736,
        148.66214192, 152.39422717, 156.16429539, 159.97319812,
        163.81871906, 167.68296294, 171.58126978, 175.51432032,
        179.4853176 , 183.47946829, 187.50021625, 191.55521759,
        195.64650875, 199.76805432, 203.91245975, 208.08280489]])

INNO2,S2  = np.meshgrid(INNO_pos,Spgrid)

splVp = interp2d(Spgrid, INNO_pos, values,  'cubic' )

innop_prob = array([2.00643263e-07, 3.71183088e-05, 2.57503446e-03, 7.93954096e-02,
       9.17992237e-01])


Sp3 = array([[  4.72608005,   4.72608005,   4.72608005,   4.72608005,
          4.72608005],
       [  4.95158936,   4.95158936,   4.95158936,   4.95158936,
          4.95158936],
       [  5.09633368,   5.09633368,   5.09633368,   5.09633368,
          5.09633368],
       [  5.17697075,   5.17697075,   5.17697075,   5.17697075,
          5.17697075],
       [  5.21225538,   5.21225538,   5.21225538,   5.21225538,
          5.21225538],
       [  5.2209424 ,   5.2209424 ,   5.2209424 ,   5.2209424 ,
          5.2209424 ],
       [  5.22166639,   5.22166639,   5.22166639,   5.22166639,
          5.22166639],
       [  5.23112659,   5.23112659,   5.23112659,   5.23112659,
          5.23112659],
       [  5.2644595 ,   5.2644595 ,   5.2644595 ,   5.2644595 ,
          5.2644595 ],
       [  5.33675701,   5.33675701,   5.33675701,   5.33675701,
          5.33675701],
       [  5.463111  ,   5.463111  ,   5.463111  ,   5.463111  ,
          5.463111  ],
       [  5.65861333,   5.65861333,   5.65861333,   5.65861333,
          5.65861333],
       [  5.93835588,   5.93835588,   5.93835588,   5.93835588,
          5.93835588],
       [  6.31803116,   6.31803116,   6.31803116,   6.31803116,
          6.31803116],
       [  6.81599595,   6.81599595,   6.81599595,   6.81599595,
          6.81599595],
       [  7.44464869,   7.44464869,   7.44464869,   7.44464869,
          7.44464869],
       [  8.21923533,   8.21923533,   8.21923533,   8.21923533,
          8.21923533],
       [  9.15500931,   9.15500931,   9.15500931,   9.15500931,
          9.15500931],
       [ 10.2671854 ,  10.2671854 ,  10.2671854 ,  10.2671854 ,
         10.2671854 ],
       [ 11.57102447,  11.57102447,  11.57102447,  11.57102447,
         11.57102447],
       [ 13.08175248,  13.08175248,  13.08175248,  13.08175248,
         13.08175248],
       [ 14.81459893,  14.81459893,  14.81459893,  14.81459893,
         14.81459893],
       [  6.21970222,   6.21970222,   6.21970222,   6.21970222,
          6.21970222],
       [  6.93325909,   6.93325909,   6.93325909,   6.93325909,
          6.93325909],
       [  7.7327826 ,   7.7327826 ,   7.7327826 ,   7.7327826 ,
          7.7327826 ],
       [  8.62316226,   8.62316226,   8.62316226,   8.62316226,
          8.62316226],
       [  9.60928761,   9.60928761,   9.60928761,   9.60928761,
          9.60928761],
       [ 10.69604815,  10.69604815,  10.69604815,  10.69604815,
         10.69604815],
       [ 11.88833341,  11.88833341,  11.88833341,  11.88833341,
         11.88833341],
       [ 13.19103289,  13.19103289,  13.19103289,  13.19103289,
         13.19103289],
       [ 14.60903612,  14.60903612,  14.60903612,  14.60903612,
         14.60903612],
       [ 16.14723261,  16.14723261,  16.14723261,  16.14723261,
         16.14723261],
       [ 17.81051188,  17.81051188,  17.81051188,  17.81051188,
         17.81051188],
       [ 19.60376345,  19.60376345,  19.60376345,  19.60376345,
         19.60376345],
       [ 21.53187684,  21.53187684,  21.53187684,  21.53187684,
         21.53187684],
       [ 23.59974156,  23.59974156,  23.59974156,  23.59974156,
         23.59974156],
       [ 25.81224713,  25.81224713,  25.81224713,  25.81224713,
         25.81224713],
       [ 28.17428307,  28.17428307,  28.17428307,  28.17428307,
         28.17428307],
       [ 30.69073889,  30.69073889,  30.69073889,  30.69073889,
         30.69073889],
       [ 33.36650412,  33.36650412,  33.36650412,  33.36650412,
         33.36650412],
       [ 36.20646827,  36.20646827,  36.20646827,  36.20646827,
         36.20646827],
       [ 39.21552085,  39.21552085,  39.21552085,  39.21552085,
         39.21552085],
       [ 42.39855139,  42.39855139,  42.39855139,  42.39855139,
         42.39855139],
       [ 45.7604494 ,  45.7604494 ,  45.7604494 ,  45.7604494 ,
         45.7604494 ],
       [ 49.30610439,  49.30610439,  49.30610439,  49.30610439,
         49.30610439],
       [ 53.0404059 ,  53.0404059 ,  53.0404059 ,  53.0404059 ,
         53.0404059 ],
       [ 56.96824343,  56.96824343,  56.96824343,  56.96824343,
         56.96824343],
       [ 61.0945065 ,  61.0945065 ,  61.0945065 ,  61.0945065 ,
         61.0945065 ],
       [ 65.42408463,  65.42408463,  65.42408463,  65.42408463,
         65.42408463],
       [ 69.96186733,  69.96186733,  69.96186733,  69.96186733,
         69.96186733],
       [ 74.71274413,  74.71274413,  74.71274413,  74.71274413,
         74.71274413],
       [ 79.68160454,  79.68160454,  79.68160454,  79.68160454,
         79.68160454],
       [ 84.87333807,  84.87333807,  84.87333807,  84.87333807,
         84.87333807],
       [ 90.29283425,  90.29283425,  90.29283425,  90.29283425,
         90.29283425],
       [ 95.9449826 ,  95.9449826 ,  95.9449826 ,  95.9449826 ,
         95.9449826 ],
       [101.83467262, 101.83467262, 101.83467262, 101.83467262,
        101.83467262],
       [107.96679384, 107.96679384, 107.96679384, 107.96679384,
        107.96679384],
       [114.34623577, 114.34623577, 114.34623577, 114.34623577,
        114.34623577],
       [120.97788794, 120.97788794, 120.97788794, 120.97788794,
        120.97788794],
       [127.86663986, 127.86663986, 127.86663986, 127.86663986,
        127.86663986],
       [135.01738104, 135.01738104, 135.01738104, 135.01738104,
        135.01738104],
       [142.43500101, 142.43500101, 142.43500101, 142.43500101,
        142.43500101],
       [150.12438927, 150.12438927, 150.12438927, 150.12438927,
        150.12438927],
       [158.09043536, 158.09043536, 158.09043536, 158.09043536,
        158.09043536],
       [166.33802878, 166.33802878, 166.33802878, 166.33802878,
        166.33802878],
       [174.87205906, 174.87205906, 174.87205906, 174.87205906,
        174.87205906],
       [183.6974157 , 183.6974157 , 183.6974157 , 183.6974157 ,
        183.6974157 ],
       [192.81898823, 192.81898823, 192.81898823, 192.81898823,
        192.81898823],
       [202.24166617, 202.24166617, 202.24166617, 202.24166617,
        202.24166617],
       [211.97033903, 211.97033903, 211.97033903, 211.97033903,
        211.97033903],
       [222.00989633, 222.00989633, 222.00989633, 222.00989633,
        222.00989633],
       [232.36522759, 232.36522759, 232.36522759, 232.36522759,
        232.36522759],
       [243.04122233, 243.04122233, 243.04122233, 243.04122233,
        243.04122233],
       [254.04277005, 254.04277005, 254.04277005, 254.04277005,
        254.04277005],
       [265.37476028, 265.37476028, 265.37476028, 265.37476028,
        265.37476028],
       [277.04208254, 277.04208254, 277.04208254, 277.04208254,
        277.04208254],
       [289.04962635, 289.04962635, 289.04962635, 289.04962635,
        289.04962635],
       [301.40228121, 301.40228121, 301.40228121, 301.40228121,
        301.40228121],
       [314.10493666, 314.10493666, 314.10493666, 314.10493666,
        314.10493666],
       [327.1624822 , 327.1624822 , 327.1624822 , 327.1624822 ,
        327.1624822 ]])



vp = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)