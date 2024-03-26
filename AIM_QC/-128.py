from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control
from Matrix_128 import Matrix1, Matrix2, Matrix3, Matrix4, vector_b
import math

def AIM(eng):

    n = 128
    input = eng.allocate_qureg(n)
    if (resource_check != 1):
        Round_constant_XOR(eng, input, 0x112233445566778899aabbccddeeff, 128)

    if (resource_check != 1):
        print('mer_exp_3 & 27 input: ')
        print_state(eng, input, 32)
    state0, state1, ancilla = mer_exp_3(eng, input) # Toffoli depth 6

    if (resource_check != 1):
        print('mer_exp_3 result')
        print_state(eng, state0, 32)

        print('mer_exp_27 result')
        print_state(eng, state1, 32)

    out = []

    out= Matrix_Mul(eng, state0, state1)

    if(resource_check != 1):
        print_state(eng, out, 32)

    Round_constant_XOR(eng, out, vector_b, 128)

    out = mer_exp_5(eng, out, ancilla) # Toffoli depth 3, Total 9, T-depth 36

    # Feed back
    feedback(eng, input, out)
    if (resource_check != 1):
        print('Ciphertext')
        print_state(eng, out, 32)

def Round_constant_XOR(eng, k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]

def Squaring(eng, input):

    new_in = []

    for i in range(64):
        new_in.append(input[i])
        new_in.append(input[64 + i])

    # x[127]
    CNOT | (input[127], new_in[126])

    CNOT | (input[127], new_in[12])
    CNOT | (input[127], new_in[6])
    CNOT | (input[127], new_in[5])
    CNOT | (input[127], new_in[2])
    CNOT | (input[127], new_in[1])
    CNOT | (input[127], new_in[0])

    # x[126]
    CNOT | (input[126], new_in[124])  # 125
    CNOT | (input[126], new_in[126])

    CNOT | (input[126], new_in[10])  # 127
    CNOT | (input[126], new_in[5])
    CNOT | (input[126], new_in[4])
    CNOT | (input[126], new_in[3])

    # x[125]
    CNOT | (input[125], new_in[122]) #123
    CNOT | (input[125], new_in[124])

    CNOT | (input[125], new_in[8]) #
    CNOT | (input[125], new_in[3])
    CNOT | (input[125], new_in[2])
    CNOT | (input[125], new_in[1])

    for i in reversed(range(61)):
        CNOT | (input[64 + i], new_in[7 + 2 * i])
        CNOT | (input[64 + i], new_in[2 + 2 * i])
        CNOT | (input[64 + i], new_in[0 + 2 * i])


    return new_in

def Matrix_Mul(eng, state0, state1):

    out0 = eng.allocate_qureg(128)
    Matrix_Mul_lower(eng, state0, out0, Matrix1)

    out1 = eng.allocate_qureg(128)
    Matrix_Mul_upper(eng, out0, out1, Matrix2)

    out2 = eng.allocate_qureg(128)
    Matrix_Mul_lower(eng, state1, out2, Matrix3)

    Matrix_Mul_upper(eng, out2, out1, Matrix4)

    return out1 # out1 = state[0]+state[1]

def Matrix_Mul_lower(eng, input, out, matrix):

    for i in range(128):
        for j in range(i+1):
            if ((matrix[i] >> j) & 1):
                CNOT | (input[j], out[i])

def Matrix_Mul_upper(eng, input, out, matrix):
    for i in range(128):
        for j in range(128-i):
            if ((matrix[i] << j) & 0x80000000000000000000000000000000):
                CNOT | (input[127-j], out[i])

def feedback(eng, input, result):

    for i in range(128):
        CNOT | (input[i], result[i])

def mer_exp_5(eng, input, ancilla):

    n = 128

    t1 = eng.allocate_qureg(n)
    copy(eng, input, t1, 128)

    t1 = Squaring_temp(eng, t1)

    if(resource_check!= 1):
       print_state(eng, t1, 32)

    count = 0

    t2 = []
    t2, count, ancilla = recursive_karatsuba(eng, t1, input, n, count, ancilla)
    Reduction(eng, t2)

    t2_copy = eng.allocate_qureg(128)
    copy(eng, t2, t2_copy, 128)

    if (resource_check != 1):
        print_state(eng, t2, 32)

    t2 = Squaring_temp(eng, t2)
    t2 = Squaring_temp(eng, t2)

    if (resource_check != 1):
        print_state(eng, t2, 32)

    t3 = []
    count = 0
    t3, count, ancilla = recursive_karatsuba(eng, t2, t2_copy, n, count, ancilla)
    Reduction(eng, t3)

    if (resource_check != 1):
        print_state(eng, t3, 32)

    t3 = Squaring_temp(eng, t3)

    count = 0
    out = []
    out, count, ancilla = recursive_karatsuba(eng, t3, input, n, count, ancilla)
    Reduction(eng, out)

    return out


def mer_exp_3(eng, input):

    n = 128

    t1 = eng.allocate_qureg(n)
    copy(eng, input, t1, 128)

    t1 = Squaring_temp(eng, t1)

    if(resource_check!= 1):
       print_state(eng, t1, 32)

    ancilla  = eng.allocate_qureg(4118)
    count = 0

    t2 = []
    t2, count, ancilla = recursive_karatsuba(eng, t1, input, n, count, ancilla)
    Reduction(eng, t2)

    if (resource_check != 1):
        print_state(eng, t2, 32)

    t2 = Squaring_temp(eng, t2)

    if (resource_check != 1):
        print_state(eng, t2, 32)

    out = []
    count = 0
    out, count, ancilla = recursive_karatsuba(eng, t2, input, n, count, ancilla)
    Reduction(eng, out)

    t3 = eng.allocate_qureg(n)
    copy(eng, out, t3, 128)

    ######### 27 #############
    t3_copy = eng.allocate_qureg(n)
    copy(eng, t3, t3_copy, 128)  # t4_copy = t2

    t3 = Squaring_temp(eng, t3)
    t3 = Squaring_temp(eng, t3)
    t3 = Squaring_temp(eng, t3)
    if (resource_check != 1):
        print('Third Sqr')
        print_state(eng, t3, 32)

    t4 = []
    count = 0
    t4, count, ancilla = recursive_karatsuba(eng, t3, t3_copy, n, count, ancilla)
    Reduction(eng, t4)

    t4_copy = eng.allocate_qureg(n)
    copy(eng, t4, t4_copy, 128)

    t4 = Squaring_temp(eng, t4)
    t4 = Squaring_temp(eng, t4)
    t4 = Squaring_temp(eng, t4)

    t4 = Squaring_temp(eng, t4)
    t4 = Squaring_temp(eng, t4)
    t4 = Squaring_temp(eng, t4)

    if (resource_check != 1):
        print('Fourth Sqr')
        print_state(eng, t4, 32)
    t5 = []
    count = 0
    t5, count, ancilla = recursive_karatsuba(eng, t4, t4_copy, n, count, ancilla)
    Reduction(eng, t5)

    t5_copy = eng.allocate_qureg(n)
    copy(eng, t5, t5_copy, 128)

    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)

    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)

    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)

    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)

    if (resource_check != 1):
        print('Fifth Sqr')
        print_state(eng, t5, 32)

    t6 = []
    count = 0
    t6, count, ancilla = recursive_karatsuba(eng, t5_copy, t5, n, count, ancilla)
    Reduction(eng, t6)

    t6 = Squaring_temp(eng, t6)
    t6 = Squaring_temp(eng, t6)
    t6 = Squaring_temp(eng, t6)

    if (resource_check != 1):
        print('Sixth Sqr')
        print_state(eng, t6, 32)
    out1 = []
    count = 0
    out1, count, ancilla = recursive_karatsuba(eng, t6, t3_copy, n, count, ancilla)
    Reduction(eng, out1)

    if (resource_check != 1):
        print('Result')
        print_state(eng, out1, 32)

    return out, out1, ancilla

def mer_exp_27(eng, input):

    n = 128

    t1 = eng.allocate_qureg(n)
    copy(eng, input, t1, 128)

    t1 = Squaring_temp(eng, t1)

    print('Fisrt Sqr')
    print_state(eng, t1, 32)
    ancilla  = eng.allocate_qureg(5000)
    count = 0

    t2 = []
    t2, count, ancilla = recursive_karatsuba(eng, t1, input, n, count, ancilla)
    Reduction(eng, t2)

    t2 = Squaring_temp(eng, t2)
    print('Second Sqr')
    print_state(eng, t2, 32)

    t3 = []
    count = 0
    t3, count, ancilla = recursive_karatsuba(eng, t2, input, n, count, ancilla)
    Reduction(eng, t3)

    t3_copy = eng.allocate_qureg(n)
    copy(eng, t3, t3_copy, 128) #t4_copy = t2

    t3 = Squaring_temp(eng, t3)
    t3 = Squaring_temp(eng, t3)
    t3 = Squaring_temp(eng, t3)

    print('Third Sqr')
    print_state(eng, t3, 32)

    t4 = []
    count = 0
    t4, count, ancilla = recursive_karatsuba(eng, t3, t3_copy, n, count, ancilla)
    Reduction(eng, t4)

    t4_copy = eng.allocate_qureg(n)
    copy(eng, t4, t4_copy, 128)

    t4 = Squaring_temp(eng, t4)

    t4 = Squaring_temp(eng, t4)
    t4 = Squaring_temp(eng, t4)
    t4 = Squaring_temp(eng, t4)
    t4 = Squaring_temp(eng, t4)
    t4 = Squaring_temp(eng, t4)

    print('Fourth Sqr')
    print_state(eng, t4, 32)
    t5 = []
    count = 0
    t5, count, ancilla = recursive_karatsuba(eng, t4, t4_copy, n, count, ancilla)
    Reduction(eng, t5)

    t5_copy = eng.allocate_qureg(n)
    copy(eng, t5, t5_copy, 128)

    t5 = Squaring_temp(eng, t5)

    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)
    t5 = Squaring_temp(eng, t5)

    print('Fifth Sqr')
    print_state(eng, t5, 32)

    t6 = []
    count = 0
    t6, count, ancilla = recursive_karatsuba(eng, t5_copy, t5, n, count, ancilla)
    Reduction(eng, t6)

    t6 = Squaring_temp(eng, t6)
    t6 = Squaring_temp(eng, t6)
    t6 = Squaring_temp(eng, t6)

    print('Sixth Sqr')
    print_state(eng, t6, 32)
    out = []
    count = 0
    out, count, ancilla = recursive_karatsuba(eng, t6, t3_copy, n, count, ancilla)
    Reduction(eng, out)

    print('Result')
    print_state(eng, out, 32)

    return out

def print_state(eng, b, len): # if b is 128-bit -> len is 32

    All(Measure) | b
    print('0x', end='')
    print_hex(eng, b, len)
    print('\n')

def print_hex(eng, qubits, len):

    for i in reversed(range(len)):
        temp = 0
        temp = temp + int(qubits[4 * i + 3]) * 8
        temp = temp + int(qubits[4 * i + 2]) * 4
        temp = temp + int(qubits[4 * i + 1]) * 2
        temp = temp + int(qubits[4 * i])

        temp = hex(temp)
        y = temp.replace("0x", "")
        print(y, end='')

def copy(eng, a, b, n):
    for i in range(n):
        CNOT | (a[i], b[i])

def Squaring_temp(eng, input):

    new_in = []
    count = 0
    temp = eng.allocate_qureg(3)

    for i in range(64):
        new_in.append(input[i])
        new_in.append(input[64 + i])

    CNOT | (input[125], temp[0])
    CNOT | (input[126], temp[1])
    CNOT | (input[127], temp[2])

    for i in reversed(range(61)):
        CNOT | (input[64 + i], new_in[7 + 2 * i])
        CNOT | (input[64 + i], new_in[2 + 2 * i])
        CNOT | (input[64 + i], new_in[0 + 2 * i])

    CNOT | (temp[0], new_in[124])
    CNOT | (temp[0], new_in[122])

    CNOT | (temp[0], new_in[8])
    CNOT | (temp[0], new_in[3])
    CNOT | (temp[0], new_in[2])
    CNOT | (temp[0], new_in[1])

    CNOT | (temp[1], new_in[126])
    CNOT | (temp[1], new_in[124])

    CNOT | (temp[1], new_in[10])
    CNOT | (temp[1], new_in[5])
    CNOT | (temp[1], new_in[4])
    CNOT | (temp[1], new_in[3])

    CNOT | (temp[2], new_in[126])

    CNOT | (temp[2], new_in[12])
    CNOT | (temp[2], new_in[6])
    CNOT | (temp[2], new_in[5])
    CNOT | (temp[2], new_in[2])
    CNOT | (temp[2], new_in[1])
    CNOT | (temp[2], new_in[0])

    return new_in


def recursive_karatsuba(eng, a, b, n, count, ancilla): #n=4

    if(n==1):
        c = eng.allocate_qubit()
        Toffoli_gate(eng, a, b, c)

        return c, count, ancilla

    c_len = 3**math.log(n, 2) #9 #3
    r_low = n//2    #2 #1

    if(n%2!=0):
        r_low = r_low +1 # n=3 -> 2, n=4 -> 2

    r_a = []
    r_b = []

    # Provide rooms and prepare operands
    r_a = ancilla[count:count + r_low]

    count = count + r_low

    r_b = ancilla[count:count + r_low]

    count = count + r_low

    with Compute(eng):
        for i in range(r_low):
            CNOT | (a[i], r_a[i])
        for i in range(n//2):
            CNOT | (a[r_low + i], r_a[i])
        for i in range(r_low):
            CNOT | (b[i], r_b[i])
        for i in range(n//2):
            CNOT | (b[r_low + i], r_b[i])

    # upper-part setting
    if(r_low == 1):
        c = eng.allocate_qureg(3)
        Toffoli_gate(eng, a[0], b[0], c[0])
        Toffoli_gate(eng, a[1], b[1], c[2])
        CNOT | (c[0], c[1])
        CNOT | (c[2], c[1])
        Toffoli_gate(eng, r_a, r_b, c[1])

        Uncompute(eng)
        return c, count, ancilla

    c_a = []
    c_b = []
    c_r = []

    c_a, count, ancilla = recursive_karatsuba(eng, a[0:r_low], b[0:r_low], r_low, count, ancilla)# 2 qubits     # 0~2
    c_b, count, ancilla = recursive_karatsuba(eng, a[r_low:n], b[r_low:n], n//2, count, ancilla)#2 qubits        # 3~5
    c_r, count, ancilla = recursive_karatsuba(eng, r_a[0:r_low], r_b[0:r_low], r_low, count, ancilla) #2qubits  # 6~8

    Uncompute(eng)

    result = []
    result = combine(eng, c_a, c_b, c_r, n)

    return result, count, ancilla

def combine(eng, a, b, r, n):
    if (n % 2 != 0):
        # n = 13########
        for i in range(n):
            CNOT | (a[i], r[i])
        for i in range(n - 2):
            CNOT | (b[i], r[i])

        for i in range(n // 2):
            CNOT | (a[n // 2 + 1 + i], r[i])
        for i in range(n // 2):
            CNOT | (b[i], r[n // 2 + 1 + i])

        out = []
        for i in range(n // 2 + 1):  # (2n-1) = n//2 + 1 + n ? / 13 = 3+1+7+?
            out.append(a[i])
        for i in range(n):
            out.append(r[i])
        for i in range((2 * n - 1) - n // 2 - 1 - n):
            out.append(b[n // 2 + i])

        return out

    half_n = int(n/2) #n=4
    for i in range(n-1):
        CNOT | (a[i], r[i])
        CNOT | (b[i], r[i])
    for i in range(half_n-1):
        CNOT | (a[half_n+i], r[i])
        CNOT | (b[i], r[half_n+i])

    result = []
    for i in range(half_n):
        result.append(a[i])
    for i in range(n-1):
        result.append(r[i])
    for i in range(half_n):
        result.append(b[half_n-1+i])

    return result

def room(eng, length):

    room = eng.allocate_qureg(length)

    return room

def copy(eng, a, b, length):
    for i in range(length):
        CNOT | (a[i], b[i])

def Reduction(eng, x):

    for i in range(10):
        CNOT | (x[128 + i], x[0 + i])
        CNOT | (x[128 + i], x[1 + i])
        CNOT | (x[128 + i], x[2 + i])
        CNOT | (x[128 + i], x[7 + i])

        CNOT | (x[128 + 10 + i], x[0 + 10 + i])
        CNOT | (x[128 + 10 + i], x[1 + 10 + i])
        CNOT | (x[128 + 10 + i], x[2 + 10 + i])
        CNOT | (x[128 + 10 + i], x[7 + 10 + i])

        CNOT | (x[128 + 20 + i], x[0 + 20 + i])
        CNOT | (x[128 + 20 + i], x[1 + 20 + i])
        CNOT | (x[128 + 20 + i], x[2 + 20 + i])
        CNOT | (x[128 + 20 + i], x[7 + 20 + i])

        CNOT | (x[128 + 30 + i], x[0 + 30 + i])
        CNOT | (x[128 + 30 + i], x[1 + 30 + i])
        CNOT | (x[128 + 30 + i], x[2 + 30 + i])
        CNOT | (x[128 + 30 + i], x[7 + 30 + i])

        CNOT | (x[128 + 40 + i], x[0 + 40 + i])
        CNOT | (x[128 + 40 + i], x[1 + 40 + i])
        CNOT | (x[128 + 40 + i], x[2 + 40 + i])
        CNOT | (x[128 + 40 + i], x[7 + 40 + i])

        CNOT | (x[128 + 50 + i], x[0 + 50 + i])
        CNOT | (x[128 + 50 + i], x[1 + 50 + i])
        CNOT | (x[128 + 50 + i], x[2 + 50 + i])
        CNOT | (x[128 + 50 + i], x[7 + 50 + i])

        CNOT | (x[128 + 60 + i], x[0 + 60 + i])
        CNOT | (x[128 + 60 + i], x[1 + 60 + i])
        CNOT | (x[128 + 60 + i], x[2 + 60 + i])
        CNOT | (x[128 + 60 + i], x[7 + 60 + i])

        CNOT | (x[128 + 70 + i], x[0 + 70 + i])
        CNOT | (x[128 + 70 + i], x[1 + 70 + i])
        CNOT | (x[128 + 70 + i], x[2 + 70 + i])
        CNOT | (x[128 + 70 + i], x[7 + 70 + i])

        CNOT | (x[128 + 80 + i], x[0 + 80 + i])
        CNOT | (x[128 + 80 + i], x[1 + 80 + i])
        CNOT | (x[128 + 80 + i], x[2 + 80 + i])
        CNOT | (x[128 + 80 + i], x[7 + 80 + i])

        CNOT | (x[128 + 90 + i], x[0 + 90 + i])
        CNOT | (x[128 + 90 + i], x[1 + 90 + i])
        CNOT | (x[128 + 90 + i], x[2 + 90 + i])
        CNOT | (x[128 + 90 + i], x[7 + 90 + i])

        CNOT | (x[128 + 100 + i], x[0 + 100 + i])
        CNOT | (x[128 + 100 + i], x[1 + 100 + i])
        CNOT | (x[128 + 100 + i], x[2 + 100 + i])
        CNOT | (x[128 + 100 + i], x[7 + 100 + i])

        CNOT | (x[128 + 110 + i], x[0 + 110 + i])
        CNOT | (x[128 + 110 + i], x[1 + 110 + i])
        CNOT | (x[128 + 110 + i], x[2 + 110 + i])
        CNOT | (x[128 + 110 + i], x[7 + 110 + i])


    CNOT | (x[128 + 120], x[0 + 120])
    CNOT | (x[128 + 120], x[1 + 120])
    CNOT | (x[128 + 120], x[2 + 120])
    CNOT | (x[128 + 120], x[7 + 120])

    #[249]
    CNOT | (x[249], x[121])
    CNOT | (x[249], x[122])
    CNOT | (x[249], x[123])

    # [250]


    CNOT | (x[254], x[0])
    CNOT | (x[254], x[1])
    CNOT | (x[254], x[2])

    CNOT | (x[254], x[5])
    CNOT | (x[254], x[6])
    CNOT | (x[254], x[12])


    CNOT | (x[250], x[1])
    CNOT | (x[250], x[2])
    CNOT | (x[250], x[3])
    CNOT | (x[250], x[8])

    CNOT | (x[253], x[4])
    CNOT | (x[253], x[5])
    CNOT | (x[253], x[6])
    CNOT | (x[253], x[11])

    CNOT | (x[251], x[2])
    CNOT | (x[251], x[3])
    CNOT | (x[251], x[4])
    CNOT | (x[251], x[9])

    CNOT | (x[249], x[0])
    CNOT | (x[249], x[1])
    CNOT | (x[249], x[2])
    CNOT | (x[249], x[7])

    CNOT | (x[250], x[122])
    CNOT | (x[250], x[123])
    CNOT | (x[250], x[124])

    # [251]
    CNOT | (x[251], x[123])
    CNOT | (x[251], x[124])
    CNOT | (x[251], x[125])

    # [252]
    CNOT | (x[252], x[124])
    CNOT | (x[252], x[125])
    CNOT | (x[252], x[126])

    CNOT | (x[252], x[3])
    CNOT | (x[252], x[4])
    CNOT | (x[252], x[5])
    CNOT | (x[252], x[10])

    # [253]
    CNOT | (x[253], x[125])
    CNOT | (x[253], x[126])
    CNOT | (x[253], x[127])


    # [253]
    CNOT | (x[254], x[126])
    CNOT | (x[254], x[127])

def Toffoli_gate(eng, a, b, c):

    if(resource_check != 1):
        Toffoli | (a, b, c)
    else:
        Tdag | a
        Tdag | b
        H | c
        CNOT | (c, a)
        T | a
        CNOT | (b, c)
        CNOT | (b, a)
        T | c
        Tdag | a
        CNOT | (b, c)
        CNOT | (c, a)
        T | a
        Tdag | c
        CNOT | (b, a)
        H | c

global resource_check
print('Generate Ciphertext...')
Simulate = ClassicalSimulator()
eng = MainEngine(Simulate)
resource_check = 0
AIM(eng)
eng.flush()

print('Estimate cost...')
Resource = ResourceCounter()
eng = MainEngine(Resource)
resource_check = 1
AIM(eng)

print(Resource)
print('\n')
eng.flush()
