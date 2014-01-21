#
# Physics 1321 Final Project
# University of Pittsburgh
# 
# YOUR NAME: 
#

import numpy as np


def initialize_shell(shell, **kwargs):
    shell.interact(kwargs.copy())

def clear(plot_1, plot_2, plot_3, stop, messages, **kwargs):
    plot_1.clear()
    plot_2.clear()
    plot_3.clear()
    messages.clear()
    stop.value = False

def f(x, s):
    return eval(s)

def findE(Es, Ef, tol, N, dx, f_string):
    Energy = np.linspace(Es, Ef, 50)
    psiEnd = np.zeros(50)
    xs = np.linspace(-20, 20, N)
    for j in range(50):
        psi = np.zeros(N)
        psi[0] = 0.00
        psi[1] = dx
        for i in range(N-2):
            V = f(xs[i], f_string)
            psi[i+2] = (-psi[i] + 2.*psi[i+1] - 2.*(Energy[j]-V)*(dx**2)*psi[i+1])
        psiEnd[j] = psi[-1]
    temp = psiEnd.copy()
    temp = abs(temp)

    if abs(np.min(psiEnd)) > tol:
        for r in range(49):
            
            if np.sign(psiEnd[r])*np.sign(psiEnd[r+1]) == -1:
                if np.sign(psiEnd[r]) == -1:
                    findE(Energy[r+1], Energy[r], tol, N, dx, f_string)
                else:
                    findE(Energy[r], Energy[r+1], tol, N, dx, f_string)
    return Energy[temp.argmin()]
    
def findE_ISW(Es, Ef, tol, N, dx):
    Energy = np.linspace(Es, Ef, 50)
    psiEnd = np.zeros(50)
    for j in range(50):
        psi = np.zeros(N)
        psi[0] = 0.00
        psi[1] = dx
        for i in range(N-2):
            psi[i+2] = (-psi[i] + 2.*psi[i+1] - 2.*(Energy[j])*(dx**2)*psi[i+1])
        psiEnd[j] = psi[-1]
    temp = psiEnd.copy()
    temp = abs(temp)

    if abs(np.min(psiEnd)) > tol:
        for r in range(49):
            
            if np.sign(psiEnd[r])*np.sign(psiEnd[r+1]) == -1:
                if np.sign(psiEnd[r]) == -1:
                    findE_ISW(Energy[r+1], Energy[r], tol, N, dx)
                else:
                    findE_ISW(Energy[r], Energy[r+1], tol, N, dx)
    
    
    return Energy[temp.argmin()]
def Normalize(psi, x_min, x_max, N):
    integral = 0.0
    dx = (x_max - x_min)/N
    for i in range(N):
        integral += dx*psi[i]
            
    return 1.0/integral

def ISW_run(ISW_N, ISW_N_E, ISW_tol, stop, messages, plot_1, plot_2, plot_3, **kwargs):
    N = ISW_N.value
    tol = ISW_tol.value
    N_states = ISW_N_E.value
    xs = np.linspace(0, 1, N)
    dx = xs[-1] - xs[-2]
    colors = ['red', 'blue', 'green', 'black', 'orange', 'purple']
    psiEnd = np.zeros(N)
    Eig = np.zeros(N_states)
    isw = 1
    E = np.linspace(0, 200, N)
    for j in range(N):
        psi = np.zeros(N)
        psi[0] = 0.01
        psi[1] = dx
        for i in range(N-2):
            psi[i+2] = (-psi[i] + 2.*psi[i+1] - 2.*E[j]*(dx**2)*psi[i+1])
        psiEnd[j] = psi[-1]
    i = 0
    for r in range(N - 2):
        if np.sign(psiEnd[r])*np.sign(psiEnd[r+1]) == -1:
            if np.sign(psiEnd[r]) == -1 and i < N_states:
                Eig[i] = findE_ISW(E[r+1], E[r], tol, N, dx)
            elif np.sign(psiEnd[r]) == 1 and i < N_states:
                Eig[i] = findE_ISW(E[r], E[r+1], tol, N, dx)
            i += 1
    for i in range(N_states):
        messages.write("Energy Eigenvalue for n = %d is %g\n" % (i+1, Eig[i]))
        messages.write("This is an error of %g\n" % (np.log10(abs(Eig[i] - (i+1)**2.*np.pi**2./2.))))
    for j in range(N_states):
        psi = np.zeros(N)
        psi[0] = 0.0
        psi[1] = dx
        psi_actual = np.sqrt(2.)*np.sin((j+1)*np.pi*xs)
        for i in range(N-2):
            #Eig[j] = 3.
            psi[i+2] = (-psi[i] + 2.*psi[i+1] - 2.*Eig[j]*(dx**2)*psi[i+1])
        psi *= np.sqrt(Normalize(psi**2, 0., 1.0, 1000))
        data = np.column_stack((xs, psi))
        data2 = np.column_stack((xs, abs(psi - psi_actual)))
        plot_1.draw_lines(data, line_color = colors[j], line_width = 2)
        plot_2.draw_lines(data2, fill_color = colors[j])
        
        
def run(N_Steps, N_E, tolerance, function, 
        stop, messages, plot_1, plot_2, plot_3, **kwargs):
    f_string = function.value
    N = N_Steps.value
    colors = ['red', 'blue', 'green', 'black', 'orange', 'purple']
    x = 0.0
    V = f(x, f_string)

    xs = np.linspace(-4, 4, N)
    dx = xs[-1] - xs[-2]
    N_states = N_E.value
    psiEnd = np.zeros(N)
    Eig = np.zeros(N_states)
    E = np.linspace(0, 10, N)
    tol = tolerance.value
    for j in range(N):
        psi = np.zeros(N)
        psi[0] = 0.0
        psi[1] = dx
        for i in range(N-2):
            V = f(xs[i], f_string)
            psi[i+2] = (-psi[i] + 2.*psi[i+1] - 2.*(E[j]-V)*(dx**2)*psi[i+1])
        psiEnd[j] = psi[-1]
    i = 0
    for r in range(N - 2):
        if np.sign(psiEnd[r])*np.sign(psiEnd[r+1]) == -1:
            if np.sign(psiEnd[r]) == -1 and i < N_states:
                Eig[i] = findE(E[r+1], E[r], tol, N, dx, f_string)
            elif np.sign(psiEnd[r]) == 1 and i < N_states:
                Eig[i] = findE(E[r], E[r+1], tol, N, dx, f_string)
            i += 1
    x = np.linspace(-.5, 3, N)
    for k in range(N):
        plot_2.draw_point((x[k], f(x[k], f_string)))
    for j in range(N_states):
        psi = np.zeros(N)
        psi[0] = 0.0
        psi[1] = dx
        for i in range(N-2):
            V = f(xs[i], f_string)
            psi[i+2] = (-psi[i] + 2.*psi[i+1] - 2.*(Eig[j]-V)*(dx**2)*psi[i+1])
        psi *= np.sqrt(Normalize(psi**2, 0., 1.0, 1000))
        data = np.column_stack((xs, psi))
        plot_1.draw_lines(data, line_color = colors[j], line_width = 2)
    messages.write('\n=== RUN =================================\n')
    for i in range(N_states):
        messages.write("Energy Eigenvalue for n = %d is %g\n" % (i, Eig[i]))
    while True:
        if stop.value:
            break
    stop.value = False
    messages.write('DONE.\n')
