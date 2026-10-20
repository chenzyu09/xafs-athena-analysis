# ============================================================
#  XAFS 作图脚本（Origin 内置 Python / originpro），Origin 2021 及以上
#  本文件由 xafs_plot.py 自动生成，数据已嵌入。
#  运行（二选一）：
#   · Origin 的 Python 控制台（提示符 >>>）：
#       exec(open(r"本文件完整路径", encoding="utf-8").read())
#   · Origin 的命令窗口（Window → Command Window）：
#       run -pyf "本文件完整路径"
#  重新运行前建议 File → New，避免与上次生成的图重名。
#  输出：每个元素 1 个数据工作簿 + FT-EXAFS 图 + XANES 图；
#        PNG 和 Origin 项目文件(.opju)保存到 OUT_DIR
#  格式：边框 3 pt；下轴刻度朝外；上轴、右轴只有边框无刻度；
#        左轴无刻度；所有文字 Arial 36 pt；标注/图例位置按字号自动计算，避免遮挡
# ============================================================
import json
import originpro as op

# ---------------- 可调参数 ----------------
OUT_DIR    = r'__OUT_DIR__'
ELEMENTS   = __ELEMENTS__
FONT       = 'Arial'
FONT_SIZE  = 36          # 所有文字字号 (pt)
FRAME_W    = 3           # 边框 / 刻度线宽 (pt)
LINE_W     = 3           # 曲线线宽 (pt)
PAGE_W_IN  = 10          # 图页宽 (inch)
PAGE_H_IN  = 8           # 图页高 (inch)
LAYER_POS  = (14, 4, 82, 74)   # 图层 左, 上, 宽, 高（占图页的 %）
BAND_HALF  = 0.2         # 阴影带半宽 (Å)
CHAR_W     = 0.62        # Arial 平均字宽 / 字号（用于估算文字宽度，偏保守）
# ------------------------------------------

LAYER_W_PT = PAGE_W_IN * 72 * LAYER_POS[2] / 100
LAYER_H_PT = PAGE_H_IN * 72 * LAYER_POS[3] / 100
TXT_H_PT   = FONT_SIZE * 1.25


def lt(obj, cmd):
    try:
        obj.lt_exec(cmd)
    except Exception as e:
        print('  [格式设置跳过]', cmd, '->', e)


def plain(s):
    """去掉 Origin 转义（下标等），用于估算宽度"""
    return s.replace('\\-(', '').replace(')', '')


def style_line(p, color, width=LINE_W):
    try:
        p.color = color
    except Exception as e:
        print('  [颜色设置跳过]', e)
    try:
        p.set_cmd(f'-w {int(width * 500)}')            # 500 = 1 pt
    except Exception:
        try:
            p.lt_exec(f'set %C -w {int(width * 500)}')
        except Exception as e:
            print('  [线宽设置跳过]', e)


def setup_page(gl):
    f = f'font({FONT})'
    lt(gl, f'page.width={PAGE_W_IN}*page.resx; page.height={PAGE_H_IN}*page.resy;')
    l, t, w, h = LAYER_POS
    lt(gl, f'layer.unit=1; layer.left={l}; layer.top={t}; layer.width={w}; layer.height={h};')
    # 边框：四边都显示，线宽 FRAME_W
    lt(gl, 'layer.x.showAxes=3; layer.y.showAxes=3;')
    lt(gl, f'layer.x.thickness={FRAME_W}; layer.y.thickness={FRAME_W};')
    lt(gl, f'layer.x2.thickness={FRAME_W}; layer.y2.thickness={FRAME_W};')
    # 刻度：下轴朝外（主+次刻度 out = 10）；左轴、上轴、右轴无刻度
    lt(gl, 'layer.x.ticks=10; layer.y.ticks=0;')
    lt(gl, 'layer.x2.ticks=0; layer.y2.ticks=0;')
    # 字体
    lt(gl, f'layer.x.label.font={f}; layer.y.label.font={f};')
    lt(gl, f'layer.x.label.pt={FONT_SIZE}; layer.y.label.pt={FONT_SIZE};')
    lt(gl, f'xb.font={f}; yl.font={f}; xb.fsize={FONT_SIZE}; yl.fsize={FONT_SIZE};')
    lt(gl, 'layer.x.showLabels=1;')                     # 只在下轴显示数字
    lt(gl, 'layer.y.showLabels=0;')                     # y 为 a.u.，不显示数字


def add_text(gl, name, text, x_right, y_bottom):
    lt(gl, f'label -n {name} {text};')
    lt(gl, f'{name}.font=font({FONT}); {name}.fsize={FONT_SIZE}; {name}.background=0;')
    # 右对齐 + 以底边定位（.x/.y 为对象中心，单位为坐标轴单位）
    lt(gl, f'{name}.x={x_right}-{name}.dx/2; {name}.y={y_bottom}+{name}.dy/2;')


# ------------------------------------------------------------------
#  FT：按字号自动确定每条曲线的偏移和标注位置，保证文字不压曲线
# ------------------------------------------------------------------
def ft_layout(d):
    R = d['R']
    curves = d['curves']                          # 顺序：上 → 下（样品最后）
    ys = [[v - c['offset'] for v in c['ft']] for c in curves]   # 归一化（未偏移）
    n = len(curves)
    xr = 6.0
    pt_per_A = LAYER_W_PT / xr
    order = list(range(n - 1, -1, -1))            # 从下往上排
    yrange = 1.6 * n
    for _ in range(30):
        u = yrange / LAYER_H_PT                   # 每 pt 对应多少 y 单位
        h, g = TXT_H_PT * u, 8 * u
        offs, labs = {}, {}
        off = 0.0
        for j, i in enumerate(order):
            wA = (len(plain(curves[i]['label'])) * CHAR_W * FONT_SIZE + 12) / pt_per_A
            x0 = xr - 0.12 - wA
            span = [y for x, y in zip(R, ys[i]) if x >= x0]
            m = max(span) if span else 0.3
            if j > 0:
                prev = order[j - 1]
                prev_top = labs[prev][1] + h + g
                mn = min(span) if span else 0.0
                off = max(offs[prev] + 1.08, prev_top - mn)
            offs[i] = off
            labs[i] = (xr - 0.12, off + m + g)
        top = max(offs[order[-1]] + 1.05, labs[order[-1]][1] + h + g)
        new = top + 0.05
        if abs(new - yrange) < 1e-3:
            break
        yrange = new
    return ys, offs, labs, -0.05, top


def plot_ft(el, d):
    curves = d['curves']
    n = len(curves)
    ys, offs, labs, ylo, yhi = ft_layout(d)
    nb = len(d['bands'])
    wks = op.new_sheet('w', lname=f'{el}_FT')
    wks.cols = 1 + 2 * n + 2 * nb
    wks.from_list(0, d['R'], lname='R', units='Å')
    for i, c in enumerate(curves):
        wks.from_list(1 + i, ys[i], lname=c['name'], units='a.u.', comments='归一化')
    for i, c in enumerate(curves):
        wks.from_list(1 + n + i, [v + offs[i] for v in ys[i]], lname=c['name'], units='a.u.',
                      comments=f'归一化+偏移 {offs[i]:.2f}')
    col = 1 + 2 * n
    for j, b in enumerate(d['bands']):
        wks.from_list(col, [b, b], lname=f'band{j+1}_x', comments=f'阴影带中心 {b} Å')
        wks.from_list(col + 1, [ylo, yhi], lname=f'band{j+1}_y')
        col += 2
    try:
        wks.cols_axis('X' + 'Y' * (2 * n) + 'XY' * nb)
    except Exception as e:
        print('  [列属性设置跳过]', e)

    gp = op.new_graph(template='line', lname=f'{el}_FT')
    gl = gp[0]
    col = 1 + 2 * n
    for j in range(nb):
        color = d['band_colors'][j]
        # 阴影带 = 一条很粗的竖线，线宽 = 带宽(2*BAND_HALF Å) 换算成 pt
        pb = gl.add_plot(wks, coly=col + 1, colx=col, type='l')
        style_line(pb, color, 2 * BAND_HALF * LAYER_W_PT / 6.0)
        col += 2
    for i, c in enumerate(curves):
        p = gl.add_plot(wks, coly=1 + n + i, colx=0, type='l')
        style_line(p, c['color'])
    setup_page(gl)
    gl.set_xlim(0, 6, 2)
    gl.set_ylim(ylo, yhi)
    gl.axis('x').title = 'R (Å)'
    gl.axis('y').title = r'FT\i(k)\+(3)χ(\i(k)) (a.u.)'
    lt(gl, f'xb.font=font({FONT}); yl.font=font({FONT}); xb.fsize={FONT_SIZE}; yl.fsize={FONT_SIZE};')
    lt(gl, 'label -r legend;')
    for i, c in enumerate(curves):
        x, y = labs[i]
        add_text(gl, f'lab{i+1}', c['label'], round(x, 3), round(y, 3))
    return gp


# ------------------------------------------------------------------
#  XANES：图例放左上角，自动抬高 y 轴上限，保证图例不压曲线
# ------------------------------------------------------------------
def plot_xanes(el, d):
    order = d['xan_order']
    cv = {c['name']: c for c in d['curves']}
    x0, x1 = d['xr']
    wks = op.new_sheet('w', lname=f'{el}_XANES')
    wks.cols = 2 * len(order)
    for i, nm in enumerate(order):
        c = cv[nm]
        wks.from_list(2 * i, c['E'], lname='Energy', units='eV', comments=c['label'])
        wks.from_list(2 * i + 1, c['mu'], lname=c['label'], units='a.u.', comments=c['label'])
    try:
        wks.cols_axis('XY' * len(order))
    except Exception as e:
        print('  [列属性设置跳过]', e)

    # 图例尺寸估算 → 求 y 上限
    maxlen = max(len(plain(cv[nm]['label'])) for nm in order)
    leg_w_pt = maxlen * CHAR_W * FONT_SIZE + 90
    leg_h_pt = len(order) * TXT_H_PT + 16
    pt_per_eV = LAYER_W_PT / (x1 - x0)
    xe = x0 + leg_w_pt / pt_per_eV + 1.0
    ylo = -0.05
    m_leg = max(v for nm in order for e, v in zip(cv[nm]['E'], cv[nm]['mu']) if x0 <= e <= xe)
    m_all = max(v for nm in order for e, v in zip(cv[nm]['E'], cv[nm]['mu']) if x0 <= e <= x1)
    frac = 1 - (leg_h_pt + 14) / LAYER_H_PT - 0.03
    yhi = max(m_all * 1.06, ylo + (m_leg - ylo) / frac)

    gp = op.new_graph(template='line', lname=f'{el}_XANES')
    gl = gp[0]
    for i, nm in enumerate(order):
        p = gl.add_plot(wks, coly=2 * i + 1, colx=2 * i, type='l')
        style_line(p, cv[nm]['color'])
    setup_page(gl)
    gl.set_xlim(x0, x1, d['xstep'])
    gl.set_ylim(ylo, round(yhi, 3))
    gl.axis('x').title = 'Energy (eV)'
    gl.axis('y').title = r'Normalized \g(m)(E) (a.u.)'
    lt(gl, f'xb.font=font({FONT}); yl.font=font({FONT}); xb.fsize={FONT_SIZE}; yl.fsize={FONT_SIZE};')
    lt(gl, f'legend.font=font({FONT}); legend.fsize={FONT_SIZE}; legend.background=0;')
    lt(gl, 'legend.x=layer.x.from+legend.dx/2+(layer.x.to-layer.x.from)*0.01;')
    lt(gl, 'legend.y=layer.y.to-legend.dy/2-(layer.y.to-layer.y.from)*0.01;')
    return gp


def main():
    pages = []
    for el in ELEMENTS:
        print(f'== {el} ==')
        d = DATA[el]
        pages += [(plot_ft(el, d), f'{el}_FT'), (plot_xanes(el, d), f'{el}_XANES')]
    if OUT_DIR:
        for gp, tag in pages:
            try:
                gp.save_fig(OUT_DIR + '\\' + tag + '.png', width=3000)
            except Exception as e:
                print('  [PNG 导出跳过]', tag, e)
        try:
            op.save(OUT_DIR + '\\XAFS_figures.opju')
        except Exception as e:
            print('  [项目保存跳过]', e)
    print('完成。图和项目文件在：', OUT_DIR)


# ---------------- 数据（由 xafs_plot.py 生成，请勿修改） ----------------
DATA_JSON = r'''__DATA__'''
DATA = json.loads(DATA_JSON)

main()
