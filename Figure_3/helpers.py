def create_scale(scale_domain, scale_range):
    d_min = scale_domain[0]
    d_max = scale_domain[1]
    # d_range = d_max - d_min

    r_min = scale_range[0]
    r_max = scale_range[1]
    r_range = r_max - r_min

    def the_scale(x):
        d = x - d_min
        d = d / d_max

        r = d * r_range
        r = r + r_min

        return r

    return the_scale

def reverse(row, key, max_sites):
    m = max_sites[row['zal']]
    return m - row[key]


def choose_actual(row, flip_list):
    if str(row['zal']) in flip_list:
        return [row['newstart'], row['newend']]
    else:
        return [row['start'], row['end']]

def assign_chr(row, zal_dict):
    return zal_dict[row['zal']]

def assign_chr_dog(row, dog_dict):
    return dog_dict[row['dog']]

def sorted_overlap(x, y):
    if y[0]>x[-1]:
        return []
    else:
        return [max(x[0], y[0]), min(x[-1], y[-1])]