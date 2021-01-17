"""
You will need:

  - overlap_drawer.py (this file)
  - chromosome_config.py (configuration file for chromosomes to plot)
  - helpers.py (set of useful custom functions)

Additional packages you will need:

 - pandas

Input files:

 - CIRCOS format input file
 - BAC format input file 
 - length data file e.g. CSL_refseq_scaffold_length_N

How to run:

 - At the bottom of this script edit (or uncomment) the lines in the `if __name__ == "__main__":` block
 - set b = an instance of `NewChromosomeDrawer`
 - NewChromooseDrawer takes 2 arguments, both of which are tuples, both of which take the format ('filename', 'format')
 - The first tuple is the data to be drawn on the left hand side of each pair, the second is drawn on the right hand side
 - format can be either 'CIRCOS' or 'BAC'
 - ensure the last line of the script is `b.plot()`

 e.g.

 b = NewChromosomeDrawer(
             ('DNAzoo_rbest_bundle250kb_gap500kb_bundle3_final_with_colours.txt', 'CIRCOS'),
             ('BAC_data_total_length_reversed_June2020.txt', 'BAC'))

 b.plot()

"""



from itertools import groupby

import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Rectangle, Patch

from pprint import pprint

from chromosome_config import CHROMOSOMES, COLOUR_DICT, ZALCAL_CIRCOS_ORDER, PLOT_ORDER, FLIP_LIST, get_csl_lengths, \
    FIRST_ORDER, SECOND_ORDER
from helpers import *

import pandas as pd

plt.rcParams["hatch.linewidth"] = 6


class NewChromosomeDrawer:
    """ takes circos formatted alignment and BAC file and plots pairs of chromosomes"""

    def __init__(self, first_file, second_file):

        self.first_file = first_file[0]
        self.first_file_type = first_file[1]
        self.second_file = second_file[0]
        self.second_file_type = second_file[1]

        self.X_DIM = 1.0
        self.Y_DIM = 4.0

        self.zca_chr_data = {}
        self.region_data = {}
        self.overlap_data = {}

        self.chr_translation_dict = {x[0]: x[1] for x in CHROMOSOMES}


        
        self.length_data = get_csl_lengths()

        if self.first_file_type == "CIRCOS":
            self.extract(self.first_file, self.first_file_type, 0)
        else:
            self.extract(self.first_file, self.first_file_type)

        if self.second_file_type == "CIRCOS":
            self.extract(self.second_file, self.second_file_type, 1)
        else:
            self.extract(self.second_file, self.second_file_type)

        self.get_overlaps()

        self.get_lengths()

        self.max_chr = max([x['length'] for k, x in self.zca_chr_data.items()])

        self.scale = create_scale([0, self.max_chr], [0, self.Y_DIM])
        self.flipped_scale = create_scale([0, self.max_chr], [self.Y_DIM, 0])

        self.dog_chr_data = {}
        for n, x in enumerate(CHROMOSOMES):
            k = x[1]
            if k == "CFA1_repeat":
                self.dog_chr_data['CFA1_repeat'] = {'colour': COLOUR_DICT["CFA1"], 'hatch': '//'}

            else:
                self.dog_chr_data[k] = {'colour': COLOUR_DICT[k], 'hatch': False}

    def extract(self, file, filetype, double_circos_order=0):

        if filetype == "CIRCOS":

            zal_chr = [x for x in range(1, 18)]
            zal_chr.append('X')
            if double_circos_order == 0:
                zal_dict = dict(zip(FIRST_ORDER, zal_chr))
            else:
                zal_dict = dict(zip(SECOND_ORDER, zal_chr))

            dog_dict = {x[0]: x[1] for x in CHROMOSOMES}

            
            flip_list = FLIP_LIST

            df = pd.read_csv(file, sep=" ", header=None)


            df2 = df[[0, 3, 4, 5]]
            df2.columns = ['dog', 'zal', 'start', 'end']
            g = df2.groupby('zal')
            max_sites = {}
            for name, group in g:
                max_sites[name] = max(group['end'])

            df2['zal_chr'] = df2.apply(lambda row: assign_chr(row, zal_dict), axis=1)
            df2['dog_chr'] = df2.apply(lambda row: assign_chr_dog(row, dog_dict), axis=1)
            df2['newstart'] = df2.apply(lambda row: reverse(row, 'end', max_sites), axis=1)
            df2['newend'] = df2.apply(lambda row: reverse(row, 'start', max_sites), axis=1)

            actuals = df2.apply(lambda row: choose_actual(row, flip_list), axis=1, result_type='expand')
            df2['actstart'] = actuals[0]
            df2['actend'] = actuals[1]
            df2.sort_values(by=['zal_chr', 'actstart'], inplace=True)


            for name, group in df2.groupby('zal'):
                self.region_data[name] = [[row['dog_chr'], row['zal'], row['actstart'], row['actend']] for index, row in
                                          group.iterrows()]

        else:

            with open(file, "r") as f:
                lines = [x.strip().split() for x in f.readlines()]

            lines.pop(0)
            scaffold_ix = 0
            chr_ix = 2
            start_ix = 9
            end_ix = 10
            supp_scaff_ix = 1

            lines = sorted(lines, key=lambda x: x[scaffold_ix], reverse=True)

            for k, g in groupby(lines, lambda x: x[scaffold_ix]):
                # chr, start, end

                this_group = list(g)

                self.region_data[k] = [
                    [x[chr_ix], x[supp_scaff_ix], int(x[start_ix]), int(x[end_ix])]
                    for x in this_group
                ]

    def get_overlaps(self):

        for k, v in self.region_data.items():
            overlaps = []
            for i in range(len(v) - 1):
                dog_name = "{}|{}".format(v[i][0], v[i + 1][0])
                if k == "NW_020874458.1":
                    print("{} - {}".format(k, dog_name))
                region = sorted_overlap([v[i][2], v[i][3]], [v[i + 1][2], v[i + 1][3]])
                if region:
                    overlaps.append({'dog_name': dog_name, 'region': region})
            self.overlap_data[k] = overlaps

    def get_lengths(self):

        for k, v in self.region_data.items():

            sorted_items = sorted(v, key=lambda x: x[3], reverse=True)
            if k in self.length_data.keys():
                self.zca_chr_data[k] = {'length': int(self.length_data[k]), 'centromere': False}
            else:
                self.zca_chr_data[k] = {'length': max(sorted_items[0][2], sorted_items[0][3]), 'centromere': False}

            for key, group in groupby(sorted_items, key=lambda x: x[1]):
                if key[-1] == "p":
                    sorted_sub_items = sorted(list(group), key=lambda x: x[2], reverse=False)
                    self.zca_chr_data[k]['centromere'] = sorted_sub_items[0][3]

    def plot_chromosome(self, chromosome, subplot_id, order):

        chr_scale = self.scale

        x_start_1 = 1 / 6
        x_end_1 = x_start_1 + 0.25

        x_start_2 = x_end_1 + 1 / 6
        x_end_2 = x_start_2 + 0.25

        if order == 0:
            x_start = x_start_1
            x_end = x_end_1
        else:
            x_start = x_start_2
            x_end = x_end_2

        this_ax = self.ax.reshape(-1)[subplot_id]
        scale_adjust = 0

        for i in self.region_data[chromosome]:

            if i[0] != 'gap':
                scaled_start = chr_scale(i[2])
                scaled_height = chr_scale(i[3]) - chr_scale(i[2])

                scaled_start += scale_adjust

                dog_chr = self.dog_chr_data[i[0]]
                if i[0] == "CFA1_repeat":
                    this_edgecolor = 'grey'
                else:
                    this_edgecolor = dog_chr['colour']
                r = Rectangle((x_start, scaled_start), x_end - x_start, scaled_height,
                              facecolor=dog_chr['colour'], edgecolor=this_edgecolor, hatch=dog_chr['hatch'])

                this_ax.add_patch(r)

        for i in self.overlap_data[chromosome]:

            if order == 0: # or order == 1:
                scaled_start = chr_scale(i['region'][0])
                scaled_height = chr_scale(i['region'][1]) - chr_scale(i['region'][0])

                dog1, dog2 = i['dog_name'].split('|')

                if 'gap' not in [dog1, dog2]:

                    dog1_colour = self.dog_chr_data[dog1]['colour']
                    dog2_colour = self.dog_chr_data[dog2]['colour']

                    r = Rectangle((x_start, scaled_start), x_end - x_start, scaled_height,
                                  hatch='|', facecolor=dog1_colour, edgecolor=dog2_colour)

                    this_ax.add_patch(r)

        chrm = self.zca_chr_data[chromosome]

        y_bottom = chr_scale(0)
        y_top = chr_scale(chrm['length'])

        if chrm['centromere']:
            y_centromere = chr_scale(chrm['centromere'])
        else:
            y_centromere = False

        radius = (x_end - x_start) / 2.0
        center_x = x_start + radius

        theta1 = 0
        theta2 = 180

        y_top_actual = y_top - radius
        y_bottom_actual = y_bottom + radius

        w1 = Wedge((center_x, y_top_actual), radius, theta1, theta2, width=0.00001, facecolor='white',
                   edgecolor='black')
        w2 = Wedge((center_x, y_bottom_actual), radius, theta2, theta1, width=0.00001, facecolor='white',
                   edgecolor='black')

        this_ax.add_patch(w1)
        this_ax.add_patch(w2)

        if y_centromere:
            w3 = Wedge((center_x, y_centromere - radius + scale_adjust), radius, theta1, theta2, width=0.00001,
                       facecolor='white',
                       edgecolor='black')
            w4 = Wedge((center_x, y_centromere + radius + scale_adjust), radius, theta2, theta1, width=0.00001,
                       facecolor='white',
                       edgecolor='black')

            this_ax.add_patch(w3)
            this_ax.add_patch(w4)

            this_ax.plot([x_start, x_start], [y_top_actual, y_centromere + radius], ls='-', color='black')
            this_ax.plot([x_end, x_end], [y_top_actual, y_centromere + radius], ls='-', color='black')

            this_ax.plot([x_start, x_start], [y_bottom_actual, y_centromere - radius], ls='-', color='black')
            this_ax.plot([x_end, x_end], [y_bottom_actual, y_centromere - radius], ls='-', color='black')

        else:
            this_ax.plot([x_start, x_start], [y_top_actual, y_bottom_actual], ls='-', color='black')
            this_ax.plot([x_end, x_end], [y_top_actual, y_bottom_actual], ls='-', color='black')

        if order == 1:  # add the label this time
            if subplot_id <= 16:
                label = str(subplot_id + 1)
            else:
                label = "X"

            this_ax.text((x_start_1 + x_end_2) / 2, y_bottom_actual - (self.Y_DIM * 0.1), label, rotation=0,
                         ha='center')


    def plot(self):

        self.fig, self.ax = plt.subplots(ncols=9, nrows=2)

        for i_ax in self.ax.reshape(-1):
            i_ax.set_xlim([0.0, self.X_DIM])
            i_ax.set_ylim([0.0, self.Y_DIM])
            i_ax.axis('off')

        subplot_id = 0

        for n, c in enumerate(PLOT_ORDER):
            if n % 2 == 0:
                order = 0

            else:
                order = 1

            self.plot_chromosome(c, subplot_id, order)

            if n % 2 != 0:
                subplot_id += 1

        legend_elements = [Patch(facecolor=v, label=k) for k, v in COLOUR_DICT.items() if k != "gap"]
        self.fig.legend(handles=legend_elements, loc='lower right', ncol=1)

        plt.show()


if __name__ == "__main__":
    # b = NewChromosomeDrawer(
    #     ('Within_ZC_colour_merged_Rbest_bundle250kb_gap500kb_bundle3_final_with_colours.txt', 'CIRCOS'),
    #     ('BAC_data_total_length_reversed.txt', 'BAC'))

    # b = NewChromosomeDrawer(
    #     ('DNAzoo_rbest_bundle250kb_gap500kb_bundle3_final_with_colours.txt', 'CIRCOS'),
    #     ('BAC_data_total_length_reversed.txt', 'BAC'))

    # b = NewChromosomeDrawer(
    #    ('Rbest_bundle250kb_gap500kb_bundle3_final_with_colours.txt', 'CIRCOS'),
    #    ('BAC_data_total_length_reversed.txt', 'BAC'))
    # b = NewChromosomeDrawer(
    #     ('bundle50kb_gap50kb_final_with_colours.txt', 'CIRCOS'),
    #     ('BAC_data_total_length_reversed.txt', 'BAC'))
    # b = NewChromosomeDrawer(
    #     ('bundle50kb_gap50kb_final_with_colours.txt', 'CIRCOS'),
    #     ('DNAzoo_rbest_bundle250kb_gap500kb_bundle3_final_with_colours.txt', 'CIRCOS'))
    
    #b = NewChromosomeDrawer(
    #         ('Rep_canfam_DNAzoo_rbest_gap500kb_bundle250kb_final_with_colours.txt', 'CIRCOS'),
    #         ('BAC_data_total_length_reversed_June2020.txt', 'BAC'))

    #b = NewChromosomeDrawer(
    #         ('canfam_ZalCalVGP_rbest_gap500kb_bundle250kb_bundle3_final_with_colours.txt', 'CIRCOS'),
    #         ('BAC_data_total_length_reversed_June2020.txt', 'BAC'))

    b = NewChromosomeDrawer(
             ('DNAzoo_rbest_bundle250kb_gap500kb_bundle3_final_with_colours.txt', 'CIRCOS'),
             ('BAC_data_total_length_reversed_June2020.txt', 'BAC'))


    # b = NewChromosomeDrawer(
    #    ('Rbest_bundle250kb_gap500kb_bundle3_final_with_colours.txt', 'CIRCOS'),
    #    ('BAC_data_total_length_reversed_June2020.txt', 'BAC'))


    # b = NewChromosomeDrawer(
    #    ("Rbest_bundle250kb_gap500kb_bundle3_final_with_colours.txt", 'CIRCOS'),
    #    ('BAC_data_total_length_reversed_June2020.txt', 'BAC'))

    b.plot()
