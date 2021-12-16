import sys, os

path = "/mnt/perkinsdata/tongue_STOmics/original/image- ssDNA and H&E/Mouse-tongue-5-Library Construction"
images = [
    "",
    path + "/Mouse-tongue-5-HE-previous section(Library Construction)_edited.tif",
    path + "/FP200000495BR_E5_edited.tif",
    path + "/Mouse-tongue-5-HE-next section(Library Construction)_edited.tif"
]
names = [
    "",
    "HnE Previous",
    "ssDNA Stain",
    "HnE After"
]

def main():
    images[0] = sys.argv[1]
    names[0] = images[0].split("/")[-1]
    prop_coords = [float(sys.argv[i]) for i in range(2,6)]

    sizes = []
    max_x = 0
    max_y = 0
    for im in images:
        stream = os.popen("identify -ping -format '%w %h ' '" + im + "'").read().split(" ")
        stream = [int(stream[0]), int(stream[1])]

        sizes.append(stream)
        if stream[1] > max_y: 
            max_x = stream[0]
            max_y = stream[1]

    print(sizes)
    print(max_x, max_y)


if __name__ == "__main__":
    main()