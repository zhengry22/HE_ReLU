import sys
import matplotlib.pyplot as plt
import numpy as np

def read_vector():
    line = sys.stdin.readline().strip()
    return list(map(float, line.split()))

def main():
    # 读取运行时间
    hint = sys.stdin.readline().strip()
    print(hint)
    run_time_line = sys.stdin.readline().strip()
    run_time = float(run_time_line.split(":")[1].strip().replace(" ms", ""))

    # 读取 x, relu 和 actual 向量
    x = read_vector()  # Read the line after "x: "
    relu = read_vector()  # Read the line after "relu: "
    actual = read_vector()  # Read the line after "actual: "

    # 计算 relu 和 actual 的绝对差值
    differences = np.abs(np.array(relu) - np.array(actual))

    # 找到最大绝对差值及其对应的索引
    max_diff_index = np.argmax(differences)
    max_diff = differences[max_diff_index]
    max_diff_x = x[max_diff_index]

    # 输出最大绝对差值及其对应的 x 值
    print(f"最大绝对差值: {max_diff:.2f}, 对应的 x 值: {max_diff_x:.2f}")

    # 绘制图像
    plt.plot(x, relu, label='ReLU')
    plt.plot(x, actual, label='Actual')

    # 生成 [-5, 5) 范围内的等距离点
    interval_points = np.arange(-5, 5, 0.5)  # 生成[-5, 5)范围内等距离的点

    # 找到这些等距离点在 x 向量中的索引并提取对应的 actual 值
    interval_y_actual = []
    for point in interval_points:
        # 找到 x 中最接近 point 的值的索引
        closest_index = np.argmin(np.abs(np.array(x) - point))
        interval_y_actual.append(actual[closest_index])

    # 在图上绘制 Actual 曲线上对应的 y 值
    plt.scatter(interval_points, interval_y_actual, color='red', label='Actual Interval Points', zorder=5)

    # 为每个点添加 (x, y) 坐标标签，往上移动数字
    offset = 0.2  # 调整数字的偏移量
    for point, y_value in zip(interval_points, interval_y_actual):
        plt.text(point, y_value + offset, f'({point:.2f}, {y_value:.2f})', ha='center', color='red', fontsize=8)

    plt.xlabel('x')
    plt.ylabel('Value')
    plt.title('Plot from C++ vectors')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
