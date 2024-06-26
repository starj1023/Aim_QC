from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control
from Matrix_256 import Matrix1, Matrix2, Matrix3, Matrix4, Matrix5, Matrix6, vector_b
import math

def AIM(eng):

    n = 256
    input = eng.allocate_qureg(n)
    if (resource_check != 1):
        Round_constant_XOR(eng, input, 0x00112233445566778899aabbccddeeff00112233445566778899aabbccddeeff, n)

    if (resource_check != 1):
        print('mer_exp_3 & 53 & 7 input: ')
        print_state(eng, input, 64)
    state0, state2, state1, ancilla = mer_exp_3_53_7(eng, input) # Toffoli depth 8

    if (resource_check != 1):
        print('mer_exp_3 result')
        print_state(eng, state0, 64)

        print('mer_exp_53 result')
        print_state(eng, state1, 64)

        print('mer_exp_7 result')
        print_state(eng, state2, 64)

    out = []
    out = Matrix_Mul(eng, state0, state1, state2)
        # print_state(eng, out, 64)

    if(resource_check != 1):
        print_state(eng, out, 32)
    #
    Round_constant_XOR(eng, out, vector_b, 256)
    out = mer_exp_5(eng, out, ancilla) # Toffoli depth 3, Total 11, T-depth 44

    if (resource_check != 1):
        print_state(eng, out, 32)
    #
    # Feed back
    feedback(eng, input, out)
    if (resource_check != 1):
        print('Ciphertext')
        print_state(eng, out, 64)

    #ff6db678 ... cda4

def Round_constant_XOR(eng, k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]

def Matrix_Mul(eng, state0, state1, state2):

    out0 = eng.allocate_qureg(256)
    Matrix_Mul_lower(eng, state0, out0, Matrix1)

    out1 = eng.allocate_qureg(256)
    Matrix_Mul_upper(eng, out0, out1, Matrix2)

    out2 = eng.allocate_qureg(256)
    Matrix_Mul_lower(eng, state1, out2, Matrix3)

    Matrix_Mul_upper(eng, out2, out1, Matrix4)

    out3 = eng.allocate_qureg(256)
    Matrix_Mul_lower(eng, state2, out3, Matrix5)

    Matrix_Mul_upper(eng, out3, out1, Matrix6)

    return out1 # out1 = state[0] + state[1] + state[2]

def Matrix_Mul_lower(eng, input, out, matrix):

    for i in range(256):
        for j in range(i+1):
            if ((matrix[i] >> j) & 1):
                CNOT | (input[j], out[i])

def Matrix_Mul_upper(eng, input, out, matrix):
    for i in range(256):
        for j in range(256-i):
            if ((matrix[i] << j) & 0x8000000000000000000000000000000000000000000000000000000000000000):
                CNOT | (input[255-j], out[i])

def feedback(eng, input, result):

    for i in range(256):
        CNOT | (input[i], result[i])

def mer_exp_3_53_7(eng, input):

    n = 256

    t1 = eng.allocate_qureg(n)
    copy(eng, input, t1, n)

    t1 = Squaring_temp(eng, t1)

    ancilla_exp53 = eng.allocate_qureg(12610)
    ancilla_exp7 = eng.allocate_qureg(12610)

    count = 0

    t2, count, ancilla_exp7 = recursive_karatsuba(eng, t1, input, n, count, ancilla_exp7)
    Reduction(eng, t2)

    t2_exp53 = eng.allocate_qureg(256) # use later
    copy(eng, t2, t2_exp53, n)

    t1_exp53 = Squaring_temp(eng, t2)

    count = 0
    t3_exp53, count, ancilla_exp7, = recursive_karatsuba(eng, t1_exp53, input, n, count, ancilla_exp7)
    Reduction(eng, t3_exp53)

    out = eng.allocate_qureg(n)
    copy(eng, t3_exp53, out, n) # exp3 is finished

    if(resource_check != 1):
        print_state(eng, out, 64)

    t3_exp53_copy = eng.allocate_qureg(n)
    copy(eng, t3_exp53, t3_exp53_copy, n)

    t1_exp53 = Squaring_temp(eng, t3_exp53)
    t1_exp53 = Squaring_temp(eng, t1_exp53)
    t1_exp53 = Squaring_temp(eng, t1_exp53)

    count = 0
    t4_exp53, count, ancilla_exp7, = recursive_karatsuba(eng, t1_exp53, t3_exp53_copy, n, count, ancilla_exp7)
    Reduction(eng, t4_exp53)

    t4_exp53_copy = eng.allocate_qureg(n)
    copy(eng, t4_exp53, t4_exp53_copy, n)

    t1_exp53 = Squaring_temp(eng, t4_exp53)
    t1_exp7 = eng.allocate_qureg(n)
    copy(eng, t1_exp53, t1_exp7, n)

    count = 0
    out1, count, ancilla_exp7, = recursive_karatsuba(eng, t1_exp7, input, n, count, ancilla_exp7)
    Reduction(eng, out1) # exp7 is finished

    if(resource_check != 1):
        print_state(eng, out1, 32)

    for i in range(5):
        t1_exp53 = Squaring_temp(eng, t1_exp53)

    count = 0
    t4_exp53, count, ancilla_exp53, = recursive_karatsuba(eng, t1_exp53, t4_exp53_copy, n, count, ancilla_exp53)
    Reduction(eng, t4_exp53)

    t4_exp53_copy = eng.allocate_qureg(n)
    copy(eng, t4_exp53, t4_exp53_copy, n)

    t1_exp53 = Squaring_temp(eng, t4_exp53)
    for i in range(11):
        t1_exp53 = Squaring_temp(eng, t1_exp53)

    count = 0
    t4_exp53, count, ancilla_exp53, = recursive_karatsuba(eng, t1_exp53, t4_exp53_copy, n, count, ancilla_exp53)
    Reduction(eng, t4_exp53)

    ##
    t4_exp53_copy = eng.allocate_qureg(n)
    copy(eng, t4_exp53, t4_exp53_copy, n)

    t1_exp53 = Squaring_temp(eng, t4_exp53)
    for i in range(23):
        t1_exp53 = Squaring_temp(eng, t1_exp53)

    count = 0
    t4_exp53, count, ancilla_exp53, = recursive_karatsuba(eng, t1_exp53, t4_exp53_copy, n, count, ancilla_exp53)
    Reduction(eng, t4_exp53)

    t1_exp53 = Squaring_temp(eng, t4_exp53)
    t1_exp53 = Squaring_temp(eng, t1_exp53)
    count = 0
    t4_exp53, count, ancilla_exp53, = recursive_karatsuba(eng, t1_exp53, t2_exp53, n, count, ancilla_exp53)
    Reduction(eng, t4_exp53)

    t1_exp53 = Squaring_temp(eng, t4_exp53)
    t1_exp53 = Squaring_temp(eng, t1_exp53)
    t1_exp53 = Squaring_temp(eng, t1_exp53)
    count = 0
    out2, count, ancilla_exp53, = recursive_karatsuba(eng, t1_exp53, t3_exp53_copy, n, count, ancilla_exp53)
    Reduction(eng, out2)

    return out, out1, out2, ancilla_exp7

def mer_exp_5(eng, input, ancilla):

    n = 256

    t1 = eng.allocate_qureg(n)
    copy(eng, input, t1, n)

    t1 = Squaring_temp(eng, t1)

    count = 0

    t2 = []
    t2, count, ancilla = recursive_karatsuba(eng, t1, input, n, count, ancilla)
    Reduction(eng, t2)

    t2_copy = eng.allocate_qureg(n)
    copy(eng, t2, t2_copy, n)

    t2 = Squaring_temp(eng, t2)
    t2 = Squaring_temp(eng, t2)

    t3 = []
    count = 0
    t3, count, ancilla = recursive_karatsuba(eng, t2, t2_copy, n, count, ancilla)
    Reduction(eng, t3)

    t3 = Squaring_temp(eng, t3)

    count = 0
    out = []
    out, count, ancilla = recursive_karatsuba(eng, t3, input, n, count, ancilla)
    Reduction(eng, out)

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
    temp = eng.allocate_qureg(5)

    new_in.append(input[0])  # 0
    new_in.append(input[254])  # 1
    new_in.append(input[1])  # 2
    new_in.append(input[255])  # 3
    new_in.append(input[2])  # 4

    for i in range(125):
        new_in.append(input[128 + i])
        new_in.append(input[3 + i])

    new_in.append(input[253])  # 4

    CNOT | (input[251], temp[0])
    CNOT | (input[252], temp[1])
    CNOT | (input[253], temp[2])
    CNOT | (input[254], temp[3])
    CNOT | (input[255], temp[4])

    for i in reversed(range(123)):
        CNOT | (input[128 + i], new_in[10 + 2 * i])
        CNOT | (input[128 + i], new_in[2 + 2 * i])
        CNOT | (input[128 + i], new_in[0 + 2 * i])
    #
    CNOT | (temp[0], new_in[10])
    CNOT | (temp[0], new_in[5])
    CNOT | (temp[0], new_in[2])
    CNOT | (temp[0], new_in[0])
    #
    CNOT | (temp[0], new_in[248])
    CNOT | (temp[0], new_in[246])
    #
    CNOT | (temp[1], new_in[12])
    CNOT | (temp[1], new_in[7])
    CNOT | (temp[1], new_in[4])
    CNOT | (temp[1], new_in[2])
    #
    CNOT | (temp[1], new_in[250])
    CNOT | (temp[1], new_in[248])

    CNOT | (temp[2], new_in[14])
    CNOT | (temp[2], new_in[9])
    CNOT | (temp[2], new_in[6])
    CNOT | (temp[2], new_in[4])
    #
    CNOT | (temp[2], new_in[252])
    CNOT | (temp[2], new_in[250])
    #

    CNOT | (temp[3], new_in[16])
    CNOT | (temp[3], new_in[8])
    CNOT | (temp[3], new_in[3])
    #
    CNOT | (temp[3], new_in[254])
    CNOT | (temp[3], new_in[252])

    CNOT | (temp[4], new_in[18])
    CNOT | (temp[4], new_in[2])
    CNOT | (temp[4], new_in[0])
    #
    CNOT | (temp[4], new_in[254])

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

def Reduction(eng, x):

    for i in range(10):
        CNOT | (x[256 + i], x[0 + i])
        CNOT | (x[256 + i], x[2 + i])
        CNOT | (x[256 + i], x[5 + i])
        CNOT | (x[256 + i], x[10 + i])

        CNOT | (x[266 + i], x[10 + i])
        CNOT | (x[266 + i], x[12 + i])
        CNOT | (x[266 + i], x[15 + i])
        CNOT | (x[266 + i], x[20 + i])

        CNOT | (x[276 + i], x[20 + i])
        CNOT | (x[276 + i], x[22 + i])
        CNOT | (x[276 + i], x[25 + i])
        CNOT | (x[276 + i], x[30 + i])

        CNOT | (x[286 + i], x[30 + i])
        CNOT | (x[286 + i], x[32 + i])
        CNOT | (x[286 + i], x[35 + i])
        CNOT | (x[286 + i], x[40 + i])

        CNOT | (x[296 + i], x[40 + i])
        CNOT | (x[296 + i], x[42 + i])
        CNOT | (x[296 + i], x[45 + i])
        CNOT | (x[296 + i], x[50 + i])

        CNOT | (x[306 + i], x[50 + i])
        CNOT | (x[306 + i], x[52 + i])
        CNOT | (x[306 + i], x[55 + i])
        CNOT | (x[306 + i], x[60 + i])

        CNOT | (x[316 + i], x[60 + i])
        CNOT | (x[316 + i], x[62 + i])
        CNOT | (x[316 + i], x[65 + i])
        CNOT | (x[316 + i], x[70 + i])

        CNOT | (x[326 + i], x[70 + i])
        CNOT | (x[326 + i], x[72 + i])
        CNOT | (x[326 + i], x[75 + i])
        CNOT | (x[326 + i], x[80 + i])

        CNOT | (x[336 + i], x[80 + i])
        CNOT | (x[336 + i], x[82 + i])
        CNOT | (x[336 + i], x[85 + i])
        CNOT | (x[336 + i], x[90 + i])

        CNOT | (x[346 + i], x[90 + i])
        CNOT | (x[346 + i], x[92 + i])
        CNOT | (x[346 + i], x[95 + i])
        CNOT | (x[346 + i], x[100 + i])

        CNOT | (x[356 + i], x[100 + i])
        CNOT | (x[356 + i], x[102 + i])
        CNOT | (x[356 + i], x[105 + i])
        CNOT | (x[356 + i], x[110 + i])

        CNOT | (x[366 + i], x[110 + i])
        CNOT | (x[366 + i], x[112 + i])
        CNOT | (x[366 + i], x[115 + i])
        CNOT | (x[366 + i], x[120 + i])

        CNOT | (x[376 + i], x[120 + i])
        CNOT | (x[376 + i], x[122 + i])
        CNOT | (x[376 + i], x[125 + i])
        CNOT | (x[376 + i], x[130 + i])

        CNOT | (x[386 + i], x[130 + i])
        CNOT | (x[386 + i], x[132 + i])
        CNOT | (x[386 + i], x[135 + i])
        CNOT | (x[386 + i], x[140 + i])

        CNOT | (x[396 + i], x[140 + i])
        CNOT | (x[396 + i], x[142 + i])
        CNOT | (x[396 + i], x[145 + i])
        CNOT | (x[396 + i], x[150 + i])

        CNOT | (x[406 + i], x[150 + i])
        CNOT | (x[406 + i], x[152 + i])
        CNOT | (x[406 + i], x[155 + i])
        CNOT | (x[406 + i], x[160 + i])

        CNOT | (x[416 + i], x[160 + i])
        CNOT | (x[416 + i], x[162 + i])
        CNOT | (x[416 + i], x[165 + i])
        CNOT | (x[416 + i], x[170 + i])

        CNOT | (x[426 + i], x[170 + i])
        CNOT | (x[426 + i], x[172 + i])
        CNOT | (x[426 + i], x[175 + i])
        CNOT | (x[426 + i], x[180 + i])

        CNOT | (x[436 + i], x[180 + i])
        CNOT | (x[436 + i], x[182 + i])
        CNOT | (x[436 + i], x[185 + i])
        CNOT | (x[436 + i], x[190 + i])

        CNOT | (x[446 + i], x[190 + i])
        CNOT | (x[446 + i], x[192 + i])
        CNOT | (x[446 + i], x[195 + i])
        CNOT | (x[446 + i], x[200 + i])

        CNOT | (x[456 + i], x[200 + i])
        CNOT | (x[456 + i], x[202 + i])
        CNOT | (x[456 + i], x[205 + i])
        CNOT | (x[456 + i], x[210 + i])

        CNOT | (x[466 + i], x[210 + i])
        CNOT | (x[466 + i], x[212 + i])
        CNOT | (x[466 + i], x[215 + i])
        CNOT | (x[466 + i], x[220 + i])

        CNOT | (x[476 + i], x[220 + i])
        CNOT | (x[476 + i], x[222 + i])
        CNOT | (x[476 + i], x[225 + i])
        CNOT | (x[476 + i], x[230 + i])

        CNOT | (x[486 + i], x[230 + i])
        CNOT | (x[486 + i], x[232 + i])
        CNOT | (x[486 + i], x[235 + i])
        CNOT | (x[486 + i], x[240 + i])

    for i in range(6):
        CNOT | (x[496 + i], x[240 + i])
        CNOT | (x[496 + i], x[242 + i])
        CNOT | (x[496 + i], x[245 + i])
        CNOT | (x[496 + i], x[250 + i])

    ####
    CNOT | (x[502], x[251])
    CNOT | (x[502], x[248])
    CNOT | (x[502], x[246])

    CNOT | (x[502], x[0])
    CNOT | (x[502], x[2])
    CNOT | (x[502], x[5])
    CNOT | (x[502], x[10])

    ####
    CNOT | (x[503], x[252])
    CNOT | (x[503], x[249])
    CNOT | (x[503], x[247])

    CNOT | (x[503], x[1])
    CNOT | (x[503], x[3])
    CNOT | (x[503], x[6])
    CNOT | (x[503], x[11])

    ####
    CNOT | (x[504], x[253])
    CNOT | (x[504], x[250])
    CNOT | (x[504], x[248])

    CNOT | (x[504], x[2])
    CNOT | (x[504], x[4])
    CNOT | (x[504], x[7])
    CNOT | (x[504], x[12])

    ####
    CNOT | (x[505], x[254])
    CNOT | (x[505], x[251])
    CNOT | (x[505], x[249])

    CNOT | (x[505], x[3])
    CNOT | (x[505], x[5])
    CNOT | (x[505], x[8])
    CNOT | (x[505], x[13])

    ####
    CNOT | (x[506], x[255])
    CNOT | (x[506], x[252])
    CNOT | (x[506], x[250])

    CNOT | (x[506], x[4])
    CNOT | (x[506], x[6])
    CNOT | (x[506], x[9])
    CNOT | (x[506], x[14])

    ####
    CNOT | (x[507], x[253])
    CNOT | (x[507], x[251])

    CNOT | (x[507], x[7])
    CNOT | (x[507], x[15])
    CNOT | (x[507], x[2])
    CNOT | (x[507], x[0])

    ####
    CNOT | (x[508], x[254])
    CNOT | (x[508], x[252])

    CNOT | (x[508], x[8])
    CNOT | (x[508], x[16])
    CNOT | (x[508], x[3])
    CNOT | (x[508], x[1])

    ####
    CNOT | (x[509], x[255])
    CNOT | (x[509], x[253])

    CNOT | (x[509], x[9])
    CNOT | (x[509], x[17])
    CNOT | (x[509], x[4])
    CNOT | (x[509], x[2])

    ####
    CNOT | (x[510], x[254])
    CNOT | (x[510], x[18])
    CNOT | (x[510], x[2])
    CNOT | (x[510], x[0])
    CNOT | (x[510], x[3])

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

#0xf6db67...da4

print('\n')
print('Estimate cost...')
Resource = ResourceCounter()
eng = MainEngine(Resource)
resource_check = 1
AIM(eng)

print(Resource)
eng.flush()
