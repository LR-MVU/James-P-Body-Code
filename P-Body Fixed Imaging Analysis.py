import os
from skimage import io
import numpy as np
import pandas as pd
import csv
from natsort import natsorted


# The write_csv function takes a file_location as the file_name, an array of data, and a list of headers
# and produces a csv file
def write_csv(file_name, array, headers):
    with open(file_name, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(array)


# The file_filter function takes a list of file paths and a list of key words that may be contained in the file names.
# It then matchs each file to the first keyword in the list it contains and puts files that don't contain any 
# keywords in their own list. A list of these split files is then returned
def file_filter(file_list,type_list):
    split_files = [[] for _ in range(len(type_list)+1)]
    for file in file_list:
        for i in range(len(type_list)):
            if type_list[i] in file:
                split_files[i].append(file)
                break
        split_files[len(type_list)].append(file)

    for i in range(len(type_list)+1):
        split_files[i] = natsorted(split_files[i])
    return split_files


# The split_FISH_QUANT_file function splits every FISH_QUANT file in a list of files given so that the data is by 
# dendrite or soma rather then for the full image. A list of numpy arrays reperesenting the mRNA information for each
# dendrite or soma is returned
def split_FISH_QUANT_file (fish_quant_files, mRNA_folder):
    # Initialize an empty list to store the FISH-QUANT data split by dendrite or soma
    split_data = []
    # Initialize an empty list to store lists containing both the name and the number of dendrites
    list_num_dendrite_per_image = []
    # Loop through the FISH-Quant files
    for file in fish_quant_files:
        # Ensure no files that are not text files get processed
        if not ".txt" in file:
            continue
        # Name of image is all but last 24 charecters (need to match P-Body and Print data streams)
        file_name = file[:-24]
        # Initialize a counter that keeps track of the number of dendrites or soma seen in image
        counter = 0
        # Read in the lines of the FISH-QUANT file
        with open(os.path.join(mRNA_folder,file), 'r', encoding = 'utf-8') as f:
            lines = [line.strip().split('\t') for i, line in enumerate(f)]
        # Initialize an empty array to accumlate data
        array_acc = []
        # Set initial skip value as True because the first few lines have imaging information that need to be skipped
        skip = True
        # Initialize start_next variable as False because we are not expecting the data to start next line
        start_next = False
        for line in lines:
            # If start_next is True indicating we are expecting data but instead we get "Cell_Start" this means that there
            # were no mRNA detected so append an empty list and increment counter
            if start_next and line[0] == "CELL_START":
                split_data.append(np.array([]))
                counter += 1
            # Set start_next back to false
            start_next = False
            
            # If skip is true and the first element of the line is not "Pos_Y" we do not expect data next so skip remains true and
            # we move to the next line
            if line[0] != "Pos_Y" and skip:
                # If the first element of the line is "CELL_END" in case of dendrite file or "Nuclues_END" in case of Soma file
                # we expect to start on mRNAs soon so start_next is set to True to test condition above
                if line[0] == "CELL_END" or line[0] == "Nucleus_END":
                    start_next = True
                continue
            # If first element is "Pos_Y" and skip is true mRNA data is coming next so we set skip to False and go to next line
            elif skip:
                skip = False
                continue
            # If skip is not true this means we are reading in data
            else:
                # If the first element of the line is "SPOTS_END" then we have reach the end of our mRNAs and we append
                # our accumulator list as a numpy array to the split_data list
                if str(line[0]) == "SPOTS_END":
                    split_data.append(np.array(array_acc))
                    # Reset accumulator to empty for next dendrite or soma
                    array_acc = []
                    # We will have to skip again until next dendrite or soma
                    skip = True
                    # Increment counter because we have added another array to list
                    counter += 1
                else:
                    # If line is data, convert every elemnt in the line list to a float and append it to the end of the accumulator array
                    array_acc.append(map(float, line))
        # If the FISH-QUANT file ends with a start_next remaining true it means that a dendrite or soma had no mRNA sports
        # so we append an empty list
        if start_next:
            split_data.append(np.array([]))
            counter += 1
        # Add the dendrite name and the number of dendrites to the data tracking list
        list_num_dendrite_per_image.append([file_name,counter])
    return split_data, list_num_dendrite_per_image


# This function takes in a folder and and excel file (.xlsx) contained in that folder and returns a list of 
# numpy arrays of the sheets in the excel file
def sheet_reader (excel_folder, excel_file):
    sheets_dict = pd.read_excel(os.path.join(excel_folder,excel_file), sheet_name=None)
    sheets_list_numpy = []
    for key in sheets_dict:
        sheets_list_numpy.append(sheets_dict[key].to_numpy())
    return sheets_list_numpy

# This function takes a image repersented as a numpy array and returns the number of pixels with 255 intensity
# on the greyscale (the number of white pixels)
def num_pixels(image):
    return np.sum(image == 255)


# The function mRNA_density takes as input the path to the dendrite print, the numpy array of P-body data, and the numpy array
# of mRNA FISH-QUANT results. 
def dendrite_calculations(dendrite, pbody_np, mRNA_fish_quant):
    # Read in the dendrite (or soma) print as a numpy array
    image = io.imread(os.path.join(skel_print_folder,dendrite))
    # Determine the number of pixels in the dendrite (or soma) print
    dendrite_pixels = num_pixels(image)
    # If the dendrite (or soma) print has no pixels return [-1,-1,-1,-1,-1] to indicate an error
    if dendrite_pixels == 0:
        return [-1,-1,-1,-1,-1]
    
    num_mRNA = mRNA_fish_quant.shape[0]

    # Get the total P-body size by summing the size of every P-body, the number of P-bodies, and the number of mRNAs
    # in all the p-bodies
    pbody_pixels = 0
    num_pbody = 0
    num_mRNA_in_pbody = 0
    for i in range(pbody_np.shape[0]):
        pbody_pixels += int(pbody_np[i,1])
        num_mRNA_in_pbody += int(pbody_np[i,2])
        num_pbody += 1
    
    # Calculate information about the P-bodies and mRNAs in the dendrites
    fraction_dendrite_pbody = pbody_pixels/dendrite_pixels
    num_pbody_per_dendrite_area = num_pbody/dendrite_pixels
    density_mRNA_dendrite = num_mRNA/dendrite_pixels
    # Use -1 values as a way to signal to user that we have a divide by zero error
    if pbody_pixels == 0:
        density_mRNA_pbody = -1
    else:
        density_mRNA_pbody = num_mRNA_in_pbody/pbody_pixels
    if num_mRNA == 0:
        fraction_mRNA_pbody = -1
    else:
        fraction_mRNA_pbody = num_mRNA_in_pbody/num_mRNA

    # Return this calculated information as a list
    return [fraction_dendrite_pbody, num_pbody_per_dendrite_area, density_mRNA_dendrite, density_mRNA_pbody,fraction_mRNA_pbody]
        

# Read in data to be analysis. All file names are sorted using natsort to ensure matching between data streams

# Read in the Print and Skeleton Images
skel_print_folder = r"INPUT_PRINTS_PATH_YOURSELF"
skel_print_files = []
for file in os.listdir(skel_print_folder):
    if os.path.isfile(os.path.join(skel_print_folder, file)):
        skel_print_files += [file]
# Split files in folder into nucleus print (ignored), skeletons (ignored), soma prints, and dendrites prints
divided_files = file_filter(skel_print_files, ["_nucprint_","_print_","_skel_","_somaprint_"])
nuc_print_files = divided_files[0]
dendrite_print_files = natsorted(divided_files[1])
skel_files = divided_files[2]
soma_print_files = natsorted(divided_files[3])

# Read in the P body Result files SOMA
p_body_results_folder_SOMA = r"INPUT_SOMA_P-BODY_PATH_YOURSELF"
p_body_results_files = []
# Ignore summary data "Soma_Averages.xlsx" file
for file in os.listdir(p_body_results_folder_SOMA):
    if os.path.isfile(os.path.join(p_body_results_folder_SOMA, file)) and file != "Soma_Averages.xlsx":
        p_body_results_files += [file]
# Split files in P-Body results folder by type
masks_excel = file_filter(p_body_results_files, [".gif",".xlsx"])
# Masks are not used for analysis
pbody_masks = masks_excel[0]
# .xlsx are the excel files with the data we need for analysis
pbody_excel_results_files_SOMA = natsorted(masks_excel[1])
divided_masks = file_filter(pbody_masks, ["mRNA","740s_mask"])
mRNA_masks_SOMA = divided_masks[0]
total_masks_SONA = divided_masks[1]
split_masks_SOMA = divided_masks[2]

# Read in the P body Result files DENDRITE (same as above but for dendrites)
p_body_results_folder_DENDRITE = r"INPUT_dendrite_P-BODY_PATH_YOURSELF"
p_body_results_files = []
for file in os.listdir(p_body_results_folder_DENDRITE):
    if os.path.isfile(os.path.join(p_body_results_folder_DENDRITE, file)) and file != "Dendrite_Averages.xlsx":
        p_body_results_files += [file]
masks_excel = file_filter(p_body_results_files, [".gif",".xlsx"])
pbody_masks = masks_excel[0]
pbody_excel_results_files_DENDRITE = natsorted(masks_excel[1])
divided_masks = file_filter(pbody_masks, ["mRNA","740s_mask"])
mRNA_masks_DENDRITE = divided_masks[0]
total_masks_DENDRITE = divided_masks[1]
split_masks_DENDRITE = divided_masks[2]

# Read in FISH-QUANT mRNA files for Soma
mRNA_soma_folder = r"INPUT_FISH_QUANT_SOMA_FOLDER_PATH_YOURSELF"
mRNA_soma_files = []
for file in os.listdir(mRNA_soma_folder):
    # Don't include info files with "_FQ_" or "Mature" in data files read in
    if os.path.isfile(os.path.join(mRNA_soma_folder, file)) and not ("_FQ_" in file or "MATURE" in file):
        mRNA_soma_files += [file]
mRNA_soma_files = natsorted(mRNA_soma_files)

# Read in FISH-QUANT mRNA files for Dendrites (same as above except for dendrites)
mRNA_dendrite_folder = r"INPUT_FISH_QUANT_DENDRITE_FOLDER_PATH_YOURSELF"
mRNA_dendrite_files = []
for file in os.listdir(mRNA_dendrite_folder):
    if os.path.isfile(os.path.join(mRNA_dendrite_folder, file)) and not ("_FQ_" in file or "MATURE" in file):
        mRNA_dendrite_files += [file]
mRNA_dendrite_files = natsorted(mRNA_dendrite_files)

# Call split_FISH_QUANT_files to split up dendrites and soma
split_FISH_QUANT_dendrite, list_num_dendrites_per_image_FISH_QUANT = split_FISH_QUANT_file(mRNA_dendrite_files,mRNA_dendrite_folder)
split_FISH_QUANT_soma, list_num_soma_per_image_FISH_QUANT = split_FISH_QUANT_file(mRNA_soma_files,mRNA_soma_folder)

# Use the sheet reader function to read in the sheets repersenting the dendrites of a single image
# While doing this create a tracking list for the number of dendrites in each image
dendrite_pbody_numpy_list = []
list_num_dendrites_per_image_pbody = []
for i in range(len(pbody_excel_results_files_DENDRITE)):
    dendrite_list = sheet_reader(p_body_results_folder_DENDRITE, pbody_excel_results_files_DENDRITE[i])
    dendrite_name = pbody_excel_results_files_DENDRITE[i][:-35]
    list_num_dendrites_per_image_pbody.append([dendrite_name, len(dendrite_list)])
    for dendrite in dendrite_list: 
        dendrite_pbody_numpy_list.append(dendrite)

# Create a list for the number of dendrite in each image for the print data stream
list_num_dendrites_per_image_print = []
name = ""
index = -1
for i in range(len(dendrite_print_files)):
    # The dendrite number is the digits after the last underscore but before the file ending (.tif)
    den_num = int(dendrite_print_files[i].split("_")[-1][:-4])
    # If dendrite number is below nine it is a single digit, no extra charector  removed
    if den_num <= 9:
        buffer = 0
    # If above 9 but below 99 the dendrite number is two digits, one extra charector removed
    elif den_num <= 99:
        buffer = 1
    # If above 99 this means the dendrite number is 3 digits and 2 extra digits must be removed
    else:
        buffer = 2
    # Set name for dendrite print by splicing the file name to remove first 4 charectors and last (17 + buffer)
    # First four char are "MAX_"
    new_name = dendrite_print_files[i][4:-(17+buffer)]
    # If name is not same as previous file then the dendrite is from a different image
    if name != new_name:
        list_num_dendrites_per_image_print.append([new_name,1])
        index += 1
        name = new_name
    # Otherwise dendrite is from same image so increment the dendrite counter for this image
    else:
        list_num_dendrites_per_image_print[index][1] += 1

# Similiar to above but with Soma files
list_num_soma_per_image_pbody = []
soma_pbody_numpy_list = []
for i in range(len(pbody_excel_results_files_SOMA)):
    soma_list = sheet_reader(p_body_results_folder_SOMA, pbody_excel_results_files_SOMA[i])
    soma_name = pbody_excel_results_files_SOMA[i][:-31]
    list_num_soma_per_image_pbody.append([soma_name, len(soma_list)])
    for soma in soma_list:
        soma_pbody_numpy_list.append(soma)

# Similiar to above but with Soma files
list_num_soma_per_image_print = []
name = ""
index = -1
for i in range(len(soma_print_files)):
    soma_num = int(soma_print_files[i].split("_")[-1][:-4])
    if soma_num <= 9:
        buffer = 0
    elif soma_num <= 99:
        buffer = 1
    else:
        buffer = 2
    new_name = soma_print_files[i][4:-(21+buffer)]
    if name != new_name:
        list_num_soma_per_image_print.append([new_name,1])
        index += 1
        name = new_name
    else:
        list_num_soma_per_image_print[index][1] += 1

# All three data streams have the same exact images with the same number of dendrites then proceed with analysis
if list_num_dendrites_per_image_print == list_num_dendrites_per_image_FISH_QUANT == list_num_dendrites_per_image_pbody:
    dendrite_results_array = []
    image_index = 0
    counter = 0
    for i in range(len(dendrite_pbody_numpy_list)):
        # Call dendrite_calculations function to produce results
        dendrite_results = dendrite_calculations(dendrite_print_files[i], dendrite_pbody_numpy_list[i], split_FISH_QUANT_dendrite[i])
        if counter < list_num_dendrites_per_image_pbody[image_index][1]:
            counter += 1
        else:
            image_index += 1
            counter = 1
        # Image name is based on the image name of the data being called
        image_name = list_num_dendrites_per_image_pbody[image_index][0]
        # Append list of results to dendrite_results_array
        dendrite_results_array.append([image_name, (i+1)] + dendrite_results)

    # Use write_csv function to create .csv file with dendrite data
    file_path_dendrite = os.path.join(os.getcwd(), "Dendrite_Results.csv")
    headers_dendrite = ["Image Name", "Dendrite Number", "P-Body Area Percentage in Dendrite", "Num P-Body per Dendrite Area", "Number of mRNA per Dendrite Area", "Num mRNA per P-Body Area", "Fraction of mRNA in P-Bodies"]
    write_csv(file_path_dendrite, dendrite_results_array, headers_dendrite)
# If there is a mismatch in data print an error message as well as the data information to allow user to correct their data
else: 
    print("ERROR: Mismatch Dendrite Data")
    print("Number of Dendrites per Image")
    print("Prints:")
    print(list_num_dendrites_per_image_print)
    print("Pbody: ")
    print(list_num_dendrites_per_image_pbody)
    print("FISH-QUANT: ")
    print(list_num_dendrites_per_image_FISH_QUANT)

# Similiar to above but with soma files instead of dendrite
if list_num_soma_per_image_print == list_num_soma_per_image_FISH_QUANT == list_num_soma_per_image_pbody:
    soma_results_array = []
    image_index = 0
    counter = 0
    for i in range(len(soma_pbody_numpy_list)):
        if counter < list_num_soma_per_image_pbody[image_index][1]:
            counter += 1
        else:
            image_index += 1
            counter = 1
        image_name = list_num_soma_per_image_pbody[image_index][0]
        soma_results = dendrite_calculations(soma_print_files[i], soma_pbody_numpy_list[i], split_FISH_QUANT_soma[i])
        soma_results_array.append([image_name, (i+1)] + soma_results)

    headers_soma = ["Image Name", "Soma Number", "P-Body Area Percentage in Soma", "Num P-Body per Soma Area", "Number of mRNA per Soma Area", "Num mRNA per P-Body Area", "Fraction of mRNA in P-Bodies"]
    file_path_soma = os.path.join(os.getcwd(), "Soma_Results.csv")
    write_csv(file_path_soma, soma_results_array, headers_soma)
else:
    print("ERROR: Mismatch SOMA Data")
    print("Number of Soma per Image")
    print("Prints:")
    print(list_num_soma_per_image_print)
    print("Pbody: ")
    print(list_num_soma_per_image_pbody)
    print("FISH-QUANT: ")
    print(list_num_soma_per_image_FISH_QUANT)