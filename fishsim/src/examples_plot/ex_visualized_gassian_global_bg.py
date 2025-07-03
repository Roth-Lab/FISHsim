import utils

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    #Generate global background
    img = utils.glob_background((400, 400), 0.01, 10000)
    print(img.shape)
    plt.imshow(img,cmap="gray")
    plt.show()
