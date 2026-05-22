"""梯形成型滤波器 (IIR 形式).

完全沿用 HYS_5data/waveform_processing/trapezoidal_filter.py 的系数构造:
  - 分母 a = [na, -2 na, na, 0, ...] (Vi)
  - 分子 b = 稀疏 (在 0, 1, na, na+1, nb, nb+1, nc, nc+1 位置)
  - scipy.signal.lfilter(b, a, x)
要求输入信号的衰减时间常数 tau_input 与 CSA 衰减一致 (用于极零相消).
"""

import numpy as np
from scipy.signal import lfilter


def make_trap_coeffs(samplerate, rise_time, flop_time, tau_input):
    """构造梯形滤波器 IIR 系数 (b, a).

    Returns
    -------
    b, a : np.ndarray
        scipy.signal.lfilter 用 — y = lfilter(b, a, x).
    """
    dt = 1.0 / samplerate
    ta = rise_time
    tb = ta + flop_time
    tc = ta + tb
    na = int(round(ta / dt))
    nb = int(round(tb / dt))
    nc = int(round(tc / dt))
    p = np.exp(-dt / tau_input)

    a = np.zeros(nc + 3)
    a[0:3] = [na, -2 * na, na]

    b = np.zeros(nc + 2)
    b[0]      = 1.0
    b[1]      = -p
    b[na]    += -1.0
    b[na+1]  += p
    b[nb]    += -1.0
    b[nb+1]  += p
    b[nc]    +=  1.0
    b[nc+1]  += -p

    return b, a


def apply_trap(x, b, a):
    """对 (T,) 或 (N, T) 信号沿最后一维做梯形滤波."""
    return lfilter(b, a, x, axis=-1)
