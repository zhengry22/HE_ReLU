import sys
import matplotlib.pyplot as plt

def read_vector():
    line = sys.stdin.readline().strip()
    return list(map(float, line.split()))

def main():
    # 读取运行时间
    run_time_line = sys.stdin.readline().strip()
    run_time = float(run_time_line.split(":")[1].strip().replace(" ms", ""))

    # 读取 x, relu 和 actual 向量
    x = read_vector()  # Read the line after "x: "
    relu = read_vector()  # Read the line after "relu: "
    actual = read_vector()  # Read the line after "actual: "

    # 输出运行时间
    print(f"运行时间: {run_time} ms")
    print(x)
    print(relu)
    print(actual)

    # 绘制图像
    plt.plot(x, relu, label='ReLU')
    plt.plot(x, actual, label='Actual')
    plt.xlabel('x')
    plt.ylabel('Value')
    plt.title('Plot from C++ vectors')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
