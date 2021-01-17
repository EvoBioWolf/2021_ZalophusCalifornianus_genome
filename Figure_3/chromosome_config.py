from collections import OrderedDict

# for the values XXX_XXX_ORDER, enter a list of scaffold ids which correspond to the chromosomes in order from 1 to n then X

# California Sealion

ZALCAL_CIRCOS_ORDER = [
    "NW_020884868.1",
    "NW_020874458.1",
    "NW_020874464.1",
    "NW_020884869.1",
    "NW_020874474.1",
    "NW_020874467.1",
    "NW_020884870.1",
    "NW_020876041.1",
    "NW_020875357.1",
    "NW_020879070.1",
    "NW_020874544.1",
    "NW_020876465.1",
    "NW_020874462.1",
    "NW_020875730.1",
    "NW_020878390.1",
    "NW_020874788.1",
    "NW_020874921.1",
    "NW_020874514.1"
]

ZALCAL_BAC_ORDER = [
    "ZCA1",
    "ZCA2",
    "ZCA3",
    "ZCA4",
    "ZCA5",
    "ZCA6",
    "ZCA7",
    "ZCA8",
    "ZCA9",
    "ZCA10",
    "ZCA11",
    "ZCA12",
    "ZCA13",
    "ZCA14",
    "ZCA15",
    "ZCA16",
    "ZCA17",
    "ZCAX"
]

# Arc Gaz

ARC_GAZ_ORDER = [
    'CAAAJK010000084.1',
    'CAAAJK010000017.1',
    'CAAAJK010000037.1',
    'CAAAJK010000049.1',
    'CAAAJK010005180.1',
    'CAAAJK010000082.1',
    'CAAAJK010000001.1',
    'CAAAJK010000044.1',
    'CAAAJK010000002.1',
    'CAAAJK010000046.1',
    'CAAAJK010002441.1',
    'CAAAJK010000028.1',
    'CAAAJK010000008.1',
    'CAAAJK010004319.1',
    'CAAAJK010000033.1',
    'CAAAJK010000018.1',
    'CAAAJK010000025.1',
    'CAAAJK010000043.1'
]

ZALCAL_DNA_ZOO_ORDER = [
    'HiC_scaffold_2',
    'HiC_scaffold_17',
    'HiC_scaffold_14',
    'HiC_scaffold_1',
    'HiC_scaffold_16',
    'HiC_scaffold_3',
    'HiC_scaffold_5',
    'HiC_scaffold_4',
    'HiC_scaffold_8',
    'HiC_scaffold_9',
    'HiC_scaffold_6',
    'HiC_scaffold_7',
    'HiC_scaffold_13',
    'HiC_scaffold_11',
    'HiC_scaffold_15',
    'HiC_scaffold_10',
    'HiC_scaffold_12',
    'HiC_scaffold_18'
]

# Specify here whichorder list corresponds to the first file (FIRST_ORDER) and which to the second file (SECOND_ORDER)


#FIRST_ORDER = ZALCAL_CIRCOS_ORDER
#FIRST_ORDER = ARC_GAZ_ORDER
FIRST_ORDER = ZALCAL_DNA_ZOO_ORDER

#SECOND_ORDER = ZALCAL_DNA_ZOO_ORDER
SECOND_ORDER = ZALCAL_BAC_ORDER

#SECOND_ORDER = ZALCAL_CIRCOS_ORDER

# PLOT_ORDER uses the first and second order variables to pair off the chromosomes as plotted

PLOT_ORDER = []

for i in range(len(FIRST_ORDER)):
    PLOT_ORDER.append(FIRST_ORDER[i])
    PLOT_ORDER.append(SECOND_ORDER[i])

# If you need to flip the order any chromosomes in the final figure you can identify by id them here

# NOTE it's reccommended that you sort out the input data and leave FLIP_LIST blank if possible

FLIP_LIST = [] 

# FLIP_LIST = ['NW_020884868.1',
#              'NW_020884869.1',
#              'NW_020874474.1',
#              'NW_020876465.1',
#              'NW_020874462.1',
#              'NW_020875730.1']

# FLIP_LIST = ['HiC_scaffold_2',
#              'HiC_scaffold_1',
#              'HiC_scaffold_16',
#              'HiC_scaffold_7',
#              'HiC_scaffold_13',
#              'HiC_scaffold_11',
#              ]

# # FLIP_LIST = [
#         "CM019802.1",  # 1
#         "CM019804.1",  # 4
#         "CM019806.1",  # 5
#         # "CM019812.1",  # 10 (Not flipping Chromosome 10)
#         "CM019813.1",  # 12
#         "CM019816.1",  # 13
#         "CM019817.1",  # 14
#         #"CM019815.1",  # 15
# ]


# mapping from pinniped chromosome labels to colour identifiers (based on the dog chromosomes)

CHROMOSOMES = [
    ("chr1", "CFA1"),
    ("chr2", "CFA2"),
    ("chr3", "CFA3"),
    ("chr4", "CFA4"),
    ("chr5", "CFA5"),
    ("chr6", "CFA6"),
    ("chr7", "CFA7"),
    ("chr8", "CFA8"),
    ("chr9", "CFA9"),
    ("chr10", "CFA10"),
    ("chr11", "CFA11"),
    ("chr12", "CFA12"),
    ("chr13", "CFA13"),
    ("chr14", "CFA14"),
    ("chr15", "CFA15"),
    ("chr16", "CFA16"),
    ("chr17", "CFA17"),
    ("chr18", "CFA18"),
    ("chr19", "CFA19"),
    ("chr20", "CFA20"),
    ("chr21", "CFA21"),
    ("chr22", "CFA22"),
    ("chr23", "CFA23"),
    ("chr24", "CFA24"),
    ("chr25", "CFA25"),
    ("chr26", "CFA26"),
    ("chr27", "CFA27"),
    ("chr28", "CFA28"),
    ("chr29", "CFA29"),
    ("chr30", "CFA30"),
    ("chr31", "CFA31"),
    ("chr32", "CFA32"),
    ("chr33", "CFA33"),
    ("chr34", "CFA34"),
    ("chr35", "CFA35"),
    ("chr36", "CFA36"),
    ("chr37", "CFA37"),
    ("chr38", "CFA38"),
    ("chrx", "CFAX"),
    ("chrX", "CFAX"),
    ("na", "CFA1_repeat"),
]


# mapping of CSL chromosome ids to chromosome numbers - used in getting the lengths

SCAF_CHR_DICT = OrderedDict(
    [("NW_020884868.1", "1"),
     ("NW_020874458.1", "2"),
     ("NW_020874464.1", "3"),
     ("NW_020884869.1", "4"),
     ("NW_020874474.1", "5"),
     ("NW_020874467.1", "6"),
     ("NW_020884870.1", "7"),
     ("NW_020876041.1", "8"),
     ("NW_020875357.1", "9"),
     ("NW_020879070.1", "10"),
     ("NW_020874544.1", "11"),
     ("NW_020876465.1", "12"),
     ("NW_020874462.1", "13"),
     ("NW_020875730.1", "14"),
     ("NW_020878390.1", "15"),
     ("NW_020874788.1", "16"),
     ("NW_020874921.1", "17"),
     ("NW_020874514.1", "X")]
)

SCAF_CHR_DICT_BAC = OrderedDict(
    [("ZCA1", "1"),
     ("ZCA2", "2"),
     ("ZCA3", "3"),
     ("ZCA4", "4"),
     ("ZCA5", "5"),
     ("ZCA6", "6"),
     ("ZCA7", "7"),
     ("ZCA8", "8"),
     ("ZCA9", "9"),
     ("ZCA10", "10"),
     ("ZCA11", "11"),
     ("ZCA12", "12"),
     ("ZCA13", "13"),
     ("ZCA14", "14"),
     ("ZCA15", "15"),
     ("ZCA16", "16"),
     ("ZCA17", "17"),
     ("ZCAX", "X")]
)



# Colours for each identifies, the tuple is the RGB value as a decimal

COLOUR_DICT = {
    "CFAX": (0.011764706, 0.011764706, 0.011764706),
    "CFA1": (1, 0, 0.964705882),
    "CFA2": (0.054901961, 0.298039216, 0.631372549),
    "CFA3": (0.835294118, 1, 0),
    "CFA4": (0.619607843, 0, 0.556862745),
    "CFA5": (0.596078431, 1, 0.321568627),
    "CFA6": (1, 0.725490196, 0.058823529),
    "CFA7": (0.415686275, 0.509803922, 0.423529412),
    "CFA8": (0.760784314, 0.549019608, 0.623529412),
    "CFA9": (0.584313725, 0, 0.22745098),
    "CFA10": (1, 0.898039216, 0.007843137),
    "CFA11": (0.980392157, 0.921568627, 0.843137255),
    "CFA12": (0, 1, 0),
    "CFA13": (0.384313725, 0.054901961, 0),
    "CFA14": (0, 0.490196078, 0.709803922),
    "CFA15": (0.003921569, 0.815686275, 1),
    "CFA16": (1, 0.576470588, 0.494117647),
    "CFA17": (0.568627451, 0.815686275, 0.796078431),
    "CFA18": (0.756862745, 1, 0.756862745),
    "CFA19": (0.745098039, 0.6, 0.439215686),
    "CFA20": (0.588235294, 0.541176471, 0.909803922),
    "CFA21": (1, 0, 0),
    "CFA22": (0, 0.082352941, 0.266666667),
    "CFA23": (0, 0.560784314, 0.611764706),
    "CFA24": (0.419607843, 0.407843137, 0.509803922),
    "CFA25": (1, 0, 0.337254902),
    "CFA26": (0.003921569, 1, 0.996078431),
    "CFA27": (1, 0.007843137, 0.615686275),
    "CFA28": (0, 0.37254902, 0.223529412), 
    "CFA29": (1, 0.454901961, 0.639215686),
    "CFA30": (0, 0.682352941, 0.494117647),
    "CFA31": (0.654901961, 0.341176471, 0.250980392),
    "CFA32": (0, 0, 1),
    "CFA33": (0.996078431, 0.537254902, 0),
    "CFA34": (0.717647059, 0.717647059, 0.717647059),
    "CFA35": (0.733333333, 0.533333333, 0),
    "CFA36": (0.37254902, 0.678431373, 0.305882353),
    "CFA37": (0.741176471, 0.776470588, 1),
    "CFA38": (0.458823529, 0.266666667, 0.694117647),
}


# Function to get the lengths of the CSL chromosomes

def get_csl_lengths():
    with open('CSL_refseq_scaffold_length_N', 'r') as f:
        length_dict = {}
        for l in f.readlines():
            ls = l.strip().split('\t')
            n = ls[0].split(' ')
            if n[0] in list(SCAF_CHR_DICT.keys()):
                length_dict[n[0]] = ls[1]

    return length_dict
