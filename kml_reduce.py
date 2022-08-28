#
#

import argparse
from rdp import rdp
from pykml import parser as pykmlParser
from pykml.factory import KML_ElementMaker as KML
from lxml import etree

import numpy as np

def coord_to_str(coordElt):
    return str(coordElt[0]) + "," + str(coordElt[1]) + "," + str(coordElt[2])


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Reduce a KML file to less points using RDP algorithm")
    parser.add_argument('infile')
    parser.add_argument('outdir')
    parser.add_argument('--epsilon', type=int, default=-1)
    parser.add_argument('--verbose', type=int, default=0)
 
    args = parser.parse_args()

    base_file_name = args.infile.split("/").pop().split(".")[0]
    try:
        in_file = open(args.infile)
    except IOError as ioe:
        print("Failed to open: '", args.infile, "'")
        print(ioe)
        exit(1)
    with in_file:
        kml_doc = pykmlParser.parse(in_file)

    rdp_coords_list = []

    for pm in kml_doc.getroot().Document.Placemark:
        print("Placemark=", pm.name)
        if (hasattr(pm, "LineString")):
            if (args.verbose > 0):
                print("Processing coordinates for: ", pm.name)
            coords_list = pm.LineString.coordinates.text.split("\n")
            for coord in coords_list:
                oneCoord = coord.split(",")
                if (len(oneCoord) > 2):
                    rdp_coords_list.append([float(oneCoord[0]), float(oneCoord[1]), float(oneCoord[2])])

    in_file.close()

    if (args.epsilon < 0):
        # The RDP algorithm uses the perpendicular distance of each point to the line.  When the
        # perpendicular distance between the line and intermediate points is less than epsilon, the
        # points are removed (reduced).
        # The loop below simply computes the distance between each consecutive points to
        # attempt to determine a reasonable value of epsilon in the RDP algorithm.
        # In my experimentation, the median is the optimal choice.
        point_distances = []
        for idx in range(len(rdp_coords_list) - 1):
            start = np.array(rdp_coords_list[idx])
            end = np.array(rdp_coords_list[idx + 1])
            dist = np.linalg.norm(end - start)
            point_distances.append(dist)

        distance_array = np.array(point_distances)
        epsilons = [[np.ma.median(distance_array), "median"],
                    [np.ma.mean(distance_array), "mean"],
                    # [np.ma.max(distance_array) / 2, "Max_div_2"],
                    # [np.ma.min(distance_array)  * 2, "min_mult_2"],
                    # [np.ma.min(distance_array) * 3, "min_mult_3"],
                    [np.ma.min(distance_array) * 4, "min_mult_4"]
                   ]
    else:
        epsilons = [args.epsilon, string(args.epsilon)]

    for epsilon in epsilons:
        reduced_kml = rdp(rdp_coords_list, epsilon[0])
        print("Epsilon =", epsilon[1], ", Reduced List has ", len(reduced_kml),
               " elements from ", len(rdp_coords_list))

        reduced_coords = "\n        ".join(map(coord_to_str, reduced_kml))

        kml_out = KML.kml(
                   KML.Document(
                        KML.Placemark(
                            KML.name(base_file_name),
                            KML.LineString(
                                KML.coordinates(reduced_coords)
                            )
                        )
                    )
        )

        outfile_name =  args.outdir + "/" + base_file_name + "_e_" + epsilon[1] + ".kml"
        try:
            out_file = open(outfile_name, "w")
        except IOError as ioe:
            print("Failed to open: '", outfile_name, "'")
            print(ioe)
            exit(1)
        with out_file:
            out_file.write("<?xml version=\"1.0\"?>\n")
            out_str = etree.tostring(kml_out, encoding="unicode", pretty_print=True)
            out_file.write(out_str)
            out_file.close()

