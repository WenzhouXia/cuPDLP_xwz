import numpy as np
import matplotlib.pyplot as plt


def plot_heatmap_from_csv(filename, resolution):
    # 读取CSV文件到数组
    data = np.loadtxt(filename)
    print("data shape: ", data.shape)
    # 检查数据长度
    if len(data) != resolution**4:
        raise ValueError("数据长度不符合预期的 resolution^4")

    # 重塑数据为 resolution^2 x resolution^2 矩阵
    matrix = data.reshape((resolution**2, resolution**2))

    # 绘制热力图
    plt.figure(figsize=(10, 10))
    plt.imshow(matrix, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.title("Heatmap of the Data")
    plt.show()
    # 保存图片
    plt.savefig("heatmap_" + str(resolution) + ".png")
    print("Heatmap saved as heatmap_" + str(resolution) + ".png")


# 使用示例
# for resolution in [8, 16, 32]:
#     filename = "/home/xiawenzzou/Documents/XiaWenzhou/OptimalTransport/cuPDLP-C-new/build/y_solution_" + \
#         str(resolution) + ".txt"  # 这是您CSV文件的路径
#     plot_heatmap_from_csv(filename, resolution)
resolution = 64  # 您的resolution值，根据实际情况调整
filename = "/home/xiawenzhou/Documents/XiaWenzhou/OptimalTransport/cuPDLP-C-new/build/y_solution_" + \
    str(resolution) + ".txt"  # 这是您CSV文件的路径
plot_heatmap_from_csv(filename, resolution)
