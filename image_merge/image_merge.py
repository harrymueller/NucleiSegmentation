import sys, os, math

path = "bin1"
images = [
    "",
    path + "/Mouse-tongue-5-HE-previous section(Library Construction)_edited.tif",
    path + "/FP200000495BR_E5_edited.tif",
    path + "/Mouse-tongue-5-HE-next section(Library Construction)_edited.tif"
]
images = [
    "",
    "tifs/FP200000495BR_E5_edited.tif",
    path + "/Malat1_spatial.png",
    path + "/Neat1_spatial.png"
]
names = [
    "",
    "HnE Previous",
    "ssDNA Stain",
    "HnE After"
]

def get_size(path):
    """
        Returns an array containing the dimensions of the given path in pixels
    """
    stream = os.popen("identify -ping -format '%w %h ' '" + path + "'").read().split(" ")
    return [int(stream[0]), int(stream[1])]

def crop(path, crop_dim, trim=False):
    """
        Crops the given path to the sizes in dim, makes temp file .TEMP.
    """
    os.system("convert '{}' -crop {}x{}+{}+{} -delete 1 '{}'".format(
        path,
        *crop_dim,
        path.replace(".", ".TEMP.")
    ))
    if trim:
        os.system("convert -trim +repage '{}' '{}'".format(
            path.replace(".", ".TEMP."),
            path.replace(".", ".TEMP.")
        ))

def resize(path, new_dims):
    os.system("convert -resize {}x{} '{}' '{}'".format(
        *new_dims,
        path.replace(".", ".TEMP."),
        path.replace(".", ".TEMP.")
    ))

def append(path1, path2, sign):
    os.system("convert {}append '{}' '{}' '{}'".format(
        sign,
        path1,
        path2,
        path1
    ))

def rename(path1, path2):
    os.system("mv '{}' '{}'".format(path1, path2))

def clean(paths):
    for p in paths:
        os.system("rm '{}'".format(p.replace(".", ".TEMP.")))

def main():
    # getting primary image details
    images[0] = sys.argv[1]
    names[0] = images[0].split("/")[-1]

    # propartional coords
    crop_dim = [float(sys.argv[i]) for i in range(2,6)]

    # crop primary image
    dims = [get_size(images[0])]

    # ratio of x to y
    xy = dims[0][1] / dims[0][0]

    crop(images[0], [
        math.ceil(crop_dim[0] * dims[0][0]),
        math.ceil(xy * crop_dim[1] * dims[0][0]),
        math.ceil(crop_dim[2] * dims[0][0]),
        math.ceil(xy * crop_dim[3] * dims[0][0])
    ])
    max_x = crop_dim[0] * dims[0][0]
    largest = images[0]

    # crop images
    for im in images[1:]:
        dims.append(get_size(im))
        new_crop = [
            math.ceil(crop_dim[0] * dims[-1][0]),
            math.ceil(xy * crop_dim[1] * dims[-1][0]),
            math.ceil(crop_dim[2] * dims[-1][0]),
            math.ceil(xy * crop_dim[3] * dims[-1][0])
        ]
        crop(im, new_crop, trim = True)

        if new_crop[0] > max_x: 
            max_x = new_crop[0]
            largest = im

    # resize images to the largest
    size = [max_x, max_x * xy]
    for im in images:
        if im != largest:
            resize(im, size)
    
    append(images[0].replace(".", ".TEMP."), 
           images[1].replace(".", ".TEMP."), 
           "+")
    append(images[2].replace(".", ".TEMP."), 
           images[3].replace(".", ".TEMP."), 
           "+")
    
    append(images[0].replace(".", ".TEMP."),
           images[2].replace(".", ".TEMP."),
           "-")

    rename(images[0].replace(".", ".TEMP."), "final_{}_{}_{}_{}.png".format(*crop_dim))
    clean(images[1:])


if __name__ == "__main__":
    main()