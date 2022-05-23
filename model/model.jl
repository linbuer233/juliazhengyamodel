#  正压原始方程模式（Barotropic primitive equation model）

#        m=20：      纬向格点数
#        n=16：      经向格点数
#        d：         网格距
#        rm：        地图放大系数
#        f：         地转参数
#        w：         工作数组
#        cla，clo：  区域中心纬度和经度
#        dt：        时间步长
#        s：         平滑系数
#        ua，ub，uc：n-1，n，n+1时间层的纬向风速
#        va，vb，vc：n-1，n，n+1时间层的经向风速
#        za，zb，zc：n-1，n，n+1时间层的位势高度
#        na：        控制12小时预报的参数
#        nb：        记录时间积分步数的参数
#        nt2=72：    判别是否积分12小时，是否该做内点平滑；
#        nt4=6：     判定是否该做边界平滑；
#        nt5=36：    判定是否该做时间平滑。 
#        zo：        为了减小重力惯性外波的波速，增加差分格式的稳定性而引入的位势高度。
#        ni:         是否进行初始风场的静力初始化。
#                    ni=0为不进行初始化，使用读入的风场和高度场；
#                    ni=1为进行初始化，需要位势高度场做初值即可。

#-------------------------------------子程序部分-----------------------------------#
# 计算放大系数和地转参数
# rk:圆锥常数
# rlq:兰勃特投影映像平面上赤道到北极点的距离
# a:地球半径
# sita:标准余纬
# psx:区域中心余纬
# r:模式中心到北极的距离
# cla：区域中心纬度
# d：格点距
function cmf(d, cla, m, n)
    rk = 0.7156
    rlq = 11423370.0
    a = 6371000.0
    conv = 57.29578
    w1 = 2.0 / rk
    sita = 30.0 / conv               # conv=180/pi
    psx = (90.0 - cla) / conv

    # 计算模式中心到北极的距离r
    cel0 = a * sin(sita) / rk
    cel = (tan(psx / 2)) / (tan(sita / 2))
    r = cel0 * cel^rk

    # 确定网格坐标点在地图坐标系中的位置
    xi0 = -(m - 1) / 2
    yj0 = r / d + (n - 1) / 2

    # 求各格点至北极点的距离rl,(xj,yi)为模式各格点在地图坐标系中的位置
    xi1 = transpose(reshape(repeat(range(0, 19, step = 1), 16), 20, 16)) # 先生成一个20列的向量，复制16行，reshape，最后转置
    yj1 = reshape(repeat(range(0, 15, step = 1), 20), 16, 20)
    xi = xi0 .+ xi1
    yj = yj0 .- yj1

    rl = (.√(xi .^ 2 .+ yj .^ 2)) .* d

    # 求放大系数rm和科氏参数f
    w2 = (rl ./ rlq) .^ w1
    sinl = (1 .- w2) ./ (w2 .+ 1)
    rm = rk .* rl ./ (a .* .√(1 .- sinl .^ 2))
    f = 1.4584e-4 .* sinl

    return rm, f
end

# 计算地转风初值的子程序
# rm 地图放大系数
# f  地转参数
# d  格点距
function cgw(za, rm, f, d, m, n)
    # 偏H 项拿 ΔH 代替 ；偏x和偏y 项拿 网格距 d 代替 
    # 四个角上的点
    ## 右下
    ua[n, m] = -(rm[n, m] * 9.8) * (za[n-1, m] - za[n, m]) / (f[n, m] * d)
    va[n, m] = (rm[n, m] * 9.8) * (za[n, m-1] - za[n, m]) / (f[n, m] * d)
    ## 右上
    ua[1, m] = -(rm[1, m] * 9.8) * (za[2, m] - za[1, m]) / (f[1, m] * d)
    va[1, m] = (rm[1, m] * 9.8) * (za[1, m-1] - za[1, m]) / (f[1, m] * d)
    ## 左上
    ua[1, 1] = -(rm[1, 1] * 9.8) * (za[2, 1] - za[1, 1]) / (f[1, 1] * d)
    va[1, 1] = (rm[1, 1] * 9.8) * (za[1, 2] - za[1, 1]) / (f[1, 1] * d)
    ## 左下
    ua[n, 1] = -(rm[n, 1] * 9.8) * (za[n-1, 1] - za[n, 1]) / (f[n, 1] * d)
    va[n, 1] = (rm[n, 1] * 9.8) * (za[n, 2] - za[n, 1]) / (f[n, 1] * d)

    # 边界点
    # 左右边界
    for j = 1:m-1:m
        ua[2:n-1, j] = -(rm[2:n-1, j] .* 9.8) ./ f[2:n-1, j] .* ((za[3:n, j] .- za[1:n-2, j]) ./ 2) ./ d
        if j == 1
            va[2:n-1, j] = (rm[2:n-1, j] .* 9.8) ./ f[2:n-1, j] .* (za[2:n-1, j+1] .- za[2:n-1, j]) ./ d
        end
        if j == m
            va[2:n-1, j] = (rm[2:n-1, j] .* 9.8) ./ f[2:n-1, j] .* (za[2:n-1, j-1] .- za[2:n-1, j]) ./ d
        end
    end
    # 上下边界
    for i = 1:n-1:n
        va[i, 2:m-1] = (rm[i, 2:m-1] .* 9.8) ./ f[i, 2:m-1] .* ((za[i, 3:m] .- za[i, 1:m-2]) ./ 2) ./ d
        if i == 1
            ua[i, 2:m-1] = -(rm[i, 2:m-1] .* 9.8) ./ f[i, 2:m-1] .* (za[i+1, 2:m-1] .- za[i, 2:m-1]) ./ d
        end
        if i == n
            ua[i, 2:m-1] = -(rm[i, 2:m-1] .* 9.8) ./ f[i, 2:m-1] .* (za[i-1, 2:m-1] .- za[i, 2:m-1]) ./ d
        end
    end

    # 区域内点
    ua[2:n-1, 2:m-1] = -(rm[2:n-1, 2:m-1] .* 9.8) ./ f[2:n-1, 2:m-1] .* ((za[3:n, 2:m-1] .- za[1:n-2, 2:m-1]) ./ 2) ./ d
    va[2:n-1, 2:m-1] = (rm[2:n-1, 2:m-1] .* 9.8) ./ f[2:n-1, 2:m-1] .* ((za[2:n-1, 3:m] .- za[2:n-1, 1:m-2]) ./ 2) ./ d
    return ua, va
end

# 赋固定边界值  给未来两个时刻 ub , vb , zb ; uc , vc , zc 的边界赋上原始场的初值
function tbv(ub, vb, zb, ua, va, za, m, n)
    # 左右边界
    ub[1:n-1:n, :] = ua[1:n-1:n, :]
    vb[1:n-1:n, :] = va[1:n-1:n, :]
    zb[1:n-1:n, :] = za[1:n-1:n, :]
    # 上下边界
    ub[:, 1:m-1:m] = ua[:, 1:m-1:m]
    vb[:, 1:m-1:m] = va[:, 1:m-1:m]
    zb[:, 1:m-1:m] = za[:, 1:m-1:m]
    return ub, vb, zb
end

# 时间积分
# rm：地图放大系数
# f： 地转参数
# d： 格点距
# dt：时间步长
# zo：为了减小重力惯性外波的波速，增加差分格式的稳定性而引入的位势高度。
function ti(ua, va, za, ub, vb, zb, uc, vc, zc, rm, f, d, dt, zo, m, n)
    c = 1 / (4 * d)
    e = -c .* rm[2:n-1, 2:m-1] .* ((ub[2:n-1, 3:m] .+ ub[2:n-1, 2:m-1]) .*
                                   (ub[2:n-1, 3:m] .- ub[2:n-1, 2:m-1]) .+
                                   (ub[2:n-1, 2:m-1] .+ ub[2:n-1, 1:m-2]) .*
                                   (ub[2:n-1, 2:m-1] .- ub[2:n-1, 1:m-2]) .+
                                   (vb[1:n-2, 2:m-1] .+ vb[2:n-1, 2:m-1]) .*
                                   (ub[2:n-1, 2:m-1] .- ub[1:n-2, 2:m-1]) .+
                                   (vb[2:n-1, 2:m-1] .+ vb[3:n, 2:m-1]) .*
                                   (ub[3:n, 2:m-1] .- ub[2:n-1, 2:m-1]) .+ 19.6 .*
                                                                           (zb[2:n-1, 3:m] .- zb[2:n-1, 1:m-2])) .+
        f[2:n-1, 2:m-1] .* vb[2:n-1, 2:m-1]
    #------------------------------------------------
    uc[2:n-1, 2:m-1] = ua[2:n-1, 2:m-1] .+ (e .* dt)
    #------------------------------------------------
    g = -c .* rm[2:n-1, 2:m-1] .* ((ub[2:n-1, 3:m] .+ ub[2:n-1, 2:m-1]) .*
                                   (vb[2:n-1, 3:m] .- vb[2:n-1, 2:m-1]) .+
                                   (ub[2:n-1, 2:m-1] .+ ub[2:n-1, 1:m-2]) .*
                                   (vb[2:n-1, 2:m-1] .- vb[2:n-1, 1:m-2]) .+
                                   (vb[1:n-2, 2:m-1] .+ vb[2:n-1, 2:m-1]) .*
                                   (vb[2:n-1, 2:m-1] .- vb[1:n-2, 2:m-1]) .+
                                   (vb[2:n-1, 2:m-1] .+ vb[3:n, 2:m-1]) .*
                                   (vb[3:n, 2:m-1] .- vb[2:n-1, 2:m-1]) .+ 19.6 .*
                                                                           (zb[3:n, 2:m-1] .- zb[1:n-2, 2:m-1])) .-
        f[2:n-1, 2:m-1] .* ub[2:n-1, 2:m-1]
    #--------------------------------------------------
    vc[2:n-1, 2:m-1] = va[2:n-1, 2:m-1] .+ g .* dt
    #--------------------------------------------------
    h = -c .* rm[2:n-1, 2:m-1] .* ((ub[2:n-1, 3:m] .+ ub[2:n-1, 2:m-1]) .*
                                   (zb[2:n-1, 3:m] .- zb[2:n-1, 2:m-1]) .+
                                   (ub[2:n-1, 2:m-1] .+ ub[2:n-1, 1:m-2]) .*
                                   (zb[2:n-1, 2:m-1] .- zb[2:n-1, 1:m-2]) .+
                                   (vb[1:n-2, 2:m-1] .+ vb[2:n-1, 2:m-1]) .*
                                   (zb[2:n-1, 2:m-1] .- zb[1:n-2, 2:m-1]) .+
                                   (vb[2:n-1, 2:m-1] .+ vb[3:n, 2:m-1]) .*
                                   (zb[3:n, 2:m-1] .- zb[2:n-1, 2:m-1]) .+ 2 .*
                                                                           (zb[2:n-1, 2:m-1] .- zo) .* (ub[2:n-1, 3:m] .- ub[2:n-1, 1:m-2] .+
                                                                                                        vb[3:n, 2:m-1] .- vb[1:n-2, 2:m-1]))

    #---------------------------------------------------
    zc[2:n-1, 2:m-1] = za[2:n-1, 2:m-1] .+ h .* dt
    #---------------------------------------------------

    return uc, vc, zc
end

# 边界平滑子程序  边界九点平滑
# a : 需要平滑的数组
# s : 平滑系数
function ssbp(a, s, m, n)
    # 上下边界
    w = a[2:n-3:n-1, 2:m-1] .+ s ./ 2 .* (1 .- s) .* (a[2:n-3:n-1, 3:m] .+ a[3:n-3:n, 2:m-1] .+ a[2:n-3:n-1, 1:m-2] .+ a[1:n-3:n-2, 2:m-1] .- 4 .* a[2:n-3:n-1, 2:m-1])
    .+s .^ 2 ./ 4 .* (a[3:n-3:n, 3:m] .+ a[3:n-3:n, 1:m-2] .+ a[1:n-3:n-2, 1:m-2] .+ a[1:n-3:n-2, 3:m] .- 4 .* a[2:n-3:n-1, 2:m-1])

    a[2:n-3:n-1, 2:m-1] = w

    # 左右边界
    w = a[3:n-2, 2:m-3:m-1] .+ s ./ 2 .* (1 .- s) .* (a[3:n-2, 1:m-3:m-2] .+ a[3:n-2, 3:m-3:m] .+ a[2:n-3, 2:m-3:m-1] .+ a[4:n-1, 2:m-3:m-1] .- 4 .* a[3:n-2, 2:m-3:m-1])
    .+s .^ 2 ./ 4 .* (a[2:n-3, 1:m-3:m-2] .+ a[4:n-1, 1:m-3:m-2] .+ a[2:n-3, 3:m-3:m] .+ a[2:n-3, 3:m-3:m] .- 4 .* a[3:n-2, 2:m-3:m-1])

    a[3:n-2, 2:m-3:m-1] = w

    return a
end

# 内点平滑子程序，区域内点五点平滑
# s: 平滑系数
# l：选择平滑的方式 0 不平滑 ；1 正平滑 ；2 正逆平滑 
function ssip(a, s, m, n, l)
    if l == 0
        return a
    end
    # l=1 时进行正平滑
    w = a[2:n-1, 2:m-1] .+ s ./ 4 .* (a[2:n-1, 3:m] .+ a[3:n, 2:m-1] .+ a[2:n-1, 1:m-2] .+ a[1:n-2, 2:m-1] .- 4 .* a[2:n-1, 2:m-1])
    a[2:n-1, 2:m-1] = w
    # l=2 时进行逆平滑
    if l == 2
        w = a[2:n-1, 2:m-1] .- s ./ 4 .* (a[2:n-1, 3:m] .+ a[3:n, 2:m-1] .+ a[2:n-1, 1:m-2] .+ a[1:n-2, 2:m-1] .- 4 .* a[2:n-1, 2:m-1])
        a[2:n-1, 2:m-1] = w
    end
    return a
end

# 时间平滑子程序
function ts(ua, ub, uc, va, vb, vc, za, zb, zc, s, m, n)
    ub[2:n-1, 2:m-1] = ub[2:n-1, 2:m-1] .+ s .* (ua[2:n-1, 2:m-1] .+ uc[2:n-1, 2:m-1] .- 2 .* ub[2:n-1, 2:m-1]) ./ 2
    vb[2:n-1, 2:m-1] = vb[2:n-1, 2:m-1] .+ s .* (va[2:n-1, 2:m-1] .+ vc[2:n-1, 2:m-1] .- 2 .* vb[2:n-1, 2:m-1]) ./ 2
    zb[2:n-1, 2:m-1] = zb[2:n-1, 2:m-1] .+ s .* (za[2:n-1, 2:m-1] .+ zc[2:n-1, 2:m-1] .- 2 .* zb[2:n-1, 2:m-1]) ./ 2
    return ub, vb, zb
end

#------------------------------------------MAIN部分-------------------------------------#
@timev begin
    # 定义一些常数
    m = 20            #纬向格点数
    n = 16            #经向格点数
    d = 300000.0      #网格距
    cla = 51.0        #区域中心纬度
    clo = 118.0       #区域中心经度
    dt = 600.0        #时间步长
    zo = 2500.0       #为了减小重力惯性外波的波速，增加差分格式的稳定性而引入的位势高度。
    s = 0.5           #平滑系数
    nt2 = 72          #判别是否积分12小时，是否该做内点平滑；
    nt4 = 6           #判定是否该做边界平滑；
    nt5 = 36          #判定是否该做时间平滑。
    c1 = dt / 2.0     #积分半步
    c2 = dt * 2.0

    using Dates
    using CSV, Tables, DataFrames
    global ua, va, za, ub, vb, zb, uc, vc, zc
    tstart = now()
    println("#-----------欢迎使用正压原始方程模式----------#")
    println("\n")
    println("读入原始场......")
    ua = Array(CSV.read("Input\\ua.dat", DataFrame, header = false, delim = ' ', ignorerepeated = true, types = Float64))
    va = Array(CSV.read("Input\\va.dat", DataFrame, header = false, delim = ' ', ignorerepeated = true, types = Float64))
    za = Array(CSV.read("Input\\za.dat", DataFrame, header = false, delim = ' ', ignorerepeated = true, types = Float64))

    # 计算放大系数和地转参数,并写入数据文件中
    print("\n", "计算每个格点上的地图放大系数和地转参数，并写入对应输出文件......", "\n")
    rm, f = cmf(d, cla, m, n)  # 计算放大系数和地转参数子程序
    CSV.write("Output\\rm.dat", Tables.table(round.(rm, digits = 6)), header = false, delim = '\t', newline = '\n')
    CSV.write("Output\\f.dat", Tables.table(round.(f, digits = 6)), header = false, delim = '\t', newline = '\n')

    # 输入参数，选择是否静力初始化
    print("\n", "静力初始化选项(0表示不进行静力初始化、1表示进行静力初始化)，请输入:", "\n")
    print("注意:如果求地转风的子程序(风场初始化)未完成，则只能输入数字0。", "\n")
    while true
        ni = parse(Int, readline()) # 键盘输入,parse 把输入转化为 Int 型
        if ni == 1
            # 计算地转风初值
            print("进行静力初始化，由高度场求出风场......", "\n")
            global ua, va = cgw(za, rm, f, d, m, n)  # 地转风计算子程序
            # CSV.write("Output\\ub.dat", Tables.table(round.(ua, digits = 6)), header = false, delim = '\t', newline = '\n')
            # CSV.write("Output\\vb.dat", Tables.table(round.(va, digits = 6)), header = false, delim = '\t', newline = '\n')
            break
        elseif ni == 0
            print("不进行静力初始化使用给出的位势高度场和风场......", "\n")
            break
        else
            print("啊!!!输入了错误字符，请重新输入0或1!")
        end
    end

    CSV.write("Output\\ub1.dat", Tables.table(round.(ua, digits = 6)), header = false, delim = '\t', newline = '\n')
    CSV.write("Output\\vb1.dat", Tables.table(round.(va, digits = 6)), header = false, delim = '\t', newline = '\n')
    # 选择参数，选择平滑方式
    print("\n", "平滑选项:(0表示不进行平滑、1表示进行正平滑、2表示进行正逆平滑)，请输入:", "\n")
    while true
        global l = parse(Int, readline())
        local a = [l == 0, l == 1, l == 2]
        if any(a)
            break
        end
        print("平滑参数选择错误，请重新选择！", "\n")
    end

    # 给ub,vb,zb;uc,vc,zc 赋值创建一个16*20的零矩阵
    ub, vb, zb, uc, vc, zc = (zeros(16, 20), zeros(16, 20), zeros(16, 20), zeros(16, 20), zeros(16, 20), zeros(16, 20))

    print("固定边界条件赋值......", "\n")
    # 边值赋值子程序
    ub, vb, zb = tbv(ub, vb, zb, ua, va, za, m, n)
    uc, vc, zc = tbv(uc, vc, zc, ua, va, za, m, n)
    # 开始预报
    print("开始12小时预报......", "\n")

    for na in range(1, 2, step = 1)
        global ua, va, za, ub, vb, zb, uc, vc, zc
        local nb = 0
        # 欧拉后差积分1小时 积分步长 Δt ≈ 10 分钟
        for nn in range(1, 6, step = 1)
            # global ua, va, za, ub, vb, zb
            ub, vb, zb = ti(ua, va, za, ua, va, za, ub, vb, zb, rm, f, d, dt, zo, m, n)
            ua, va, za = ti(ua, va, za, ub, vb, zb, ua, va, za, rm, f, d, dt, zo, m, n)
            nb += 1
        end
        # 边界平滑
        za = ssbp(za, s, m, n)
        ua = ssbp(ua, s, m, n)
        va = ssbp(va, s, m, n)

        # 前差积分半步
        ub, vb, zb = ti(ua, va, za, ua, va, za, ub, vb, zb, rm, f, d, c1, zo, m, n)
        # 中央差积分半步
        uc, vc, zc = ti(ua, va, za, ub, vb, zb, uc, vc, zc, rm, f, d, dt, zo, m, n)
        nb += 1
        # 交换数组
        ub = uc
        vb = vc
        zb = zc

        # 中央差积分一步，共积分11个小时
        for nn = 1:66
            uc, vc, zc = ti(ua, va, za, ub, vb, zb, uc, vc, zc, rm, f, d, c2, zo, m, n)
            nb += 1
            #打印积分步数，na为大循环，nb为小循环
            print("大循环na=", na, "    小循环nb=", nb, "\n")
            if nb == nt2 # 判别是否积分12小时，是否该做内点平滑；
                break
            end

            # 判断是否做边界平滑
            if (floor(nb / nt4) * nt4) == nb
                zc = ssbp(zc, s, m, n)
                uc = ssbp(uc, s, m, n)
                vc = ssbp(vc, s, m, n)
            else
                # 判断是否做时间平滑
                a = [nb == nt5, nb == nt5 + 1]
                if any(a)
                    # 时间平滑
                    ub, vb, zb = ts(ua, ub, uc, va, vb, vc, za, zb, zc, s, m, n)
                else
                    # 交换数组
                    ua = ub
                    va = vb
                    za = zb
                    ub = uc
                    vb = vc
                    zb = zc
                end
            end
        end
        # print(nb)
        if nb == nt2
            # 内点平滑，五点平滑
            zc = ssip(zc, s, m, n, l)
            uc = ssip(uc, s, m, n, l)
            vc = ssip(vc, s, m, n, l)
            # 打印积分步数，na为大循环，nb为小循环
            print("大循环na=", na, "    小循环nb=", nb, "\n")
            # 交换数组
            ua = uc
            va = vc
            za = zc
        end
    end
    # 存放预报结果
    print("输出预报结果.......", "\n")
    CSV.write("Output\\uc.dat", Tables.table(round.(uc, digits = 6)), header = false, delim = '\t', newline = '\n')
    CSV.write("Output\\vc.dat", Tables.table(round.(vc, digits = 6)), header = false, delim = '\t', newline = '\n')
    CSV.write("Output\\zc.dat", Tables.table(zc), header = false, delim = '\t')
    print("！！！预报结束！！！")
    tend = now()
    print(tend - tstart)
end
