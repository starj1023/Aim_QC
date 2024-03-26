from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control
from Matrix_192 import Matrix1, Matrix2, Matrix3, Matrix4, vector_b

import math


def Round_constant_XOR(eng, k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]

def Matrix_Mul(eng, state0, state1):

    out0 = eng.allocate_qureg(192)
    Matrix_Mul_lower(eng, state0, out0, Matrix1)

    out1 = eng.allocate_qureg(192)
    Matrix_Mul_upper(eng, out0, out1, Matrix2)

    out2 = eng.allocate_qureg(192)
    Matrix_Mul_lower(eng, state1, out2, Matrix3)

    Matrix_Mul_upper(eng, out2, out1, Matrix4)

    return out1 # out1 = state[0]+state[1]

def Matrix_Mul_lower(eng, input, out, matrix):

    for i in range(192):
        for j in range(i+1):
            if ((matrix[i] >> j) & 1):
                CNOT | (input[j], out[i])

def Matrix_Mul_upper(eng, input, out, matrix):
    for i in range(192):
        for j in range(192-i):
            if ((matrix[i] << j) & 0x800000000000000000000000000000000000000000000000):
                CNOT | (input[191-j], out[i])

def feedback(eng, input, result):

    for i in range(192):
        CNOT | (input[i], result[i])

def AIM(eng):

    n = 192
    input = eng.allocate_qureg(n)
    if (resource_check != 1):
        Round_constant_XOR(eng, input, 0x00112233445566778899aabbccddeeff0011223344556677, n)

    if (resource_check != 1):
        print('mer_exp_5 & 29 input: ')
        print_state(eng, input, 48)
    state0, state1, ancilla = mer_exp_29_5(eng, input) # Toffoli depth 8

    if (resource_check != 1):
        print('mer_exp_5 result')
        print_state(eng, state0, 48)

        print('mer_exp_29 result')
        print_state(eng, state1, 48)

    out = []
    out= Matrix_Mul(eng, state0, state1)

    Round_constant_XOR(eng, out, vector_b, 192)
    out = mer_exp_7(eng, out, ancilla) # Toffoli depth 4, Total 12 -> T-depth 48

    # Feed back
    feedback(eng, input, out)
    if (resource_check != 1):
        print('Ciphertext')
        print_state(eng, out, 48)

def mer_exp_7(eng, input, ancilla):

    n = 192

    t1 = eng.allocate_qureg(n)
    copy(eng, input, t1, n)

    t1 = Squaring_temp(eng, t1)
    count = 0
    t1, count, ancilla = recursive_karatsuba(eng, t1, input, n, count, ancilla)
    Reduction(eng, t1)

    t1 = Squaring_temp(eng, t1)

    count = 0
    t2, count, ancilla = recursive_karatsuba(eng, t1, input, n, count, ancilla)
    Reduction(eng, t2)

    t2_copy = eng.allocate_qureg(n)
    copy(eng, t2, t2_copy, n)

    t1 = Squaring_temp(eng, t2)
    t1 = Squaring_temp(eng, t1)
    t1 = Squaring_temp(eng, t1)

    count = 0
    t1, count, ancilla = recursive_karatsuba(eng, t1, t2_copy, n, count, ancilla)
    Reduction(eng, t1)

    t1 = Squaring_temp(eng, t1)

    count = 0
    out = []
    out, count, ancilla = recursive_karatsuba(eng, t1, input, n, count, ancilla)
    Reduction(eng, out)

    return out

def mer_exp_29_5(eng, input):

    n = 192

    t1 = eng.allocate_qureg(n)
    copy(eng, input, t1, n)

    t1 = Squaring_temp(eng, t1)

    ancilla_exp5  = eng.allocate_qureg(9822)
    ancilla_exp29 = eng.allocate_qureg(9822)

    count = 0
    t1, count, ancilla_exp5 = recursive_karatsuba(eng, t1, input, n, count, ancilla_exp5)
    Reduction(eng, t1)

    t2_copy = eng.allocate_qureg(n)
    copy(eng, t1, t2_copy, n)


    t1_exp29 = Squaring_temp(eng, t1) # t1_exp29 = s -> m -> s

    t1_exp5 = eng.allocate_qureg(n)
    copy(eng, t1_exp29, t1_exp5, n) # out = sqr , mul, sqr

    if (resource_check != 1):
        print_state(eng, t1_exp5, 48) #sms

    #########
    t1_exp5 = Squaring_temp(eng, t1_exp5)
    count = 0
    t1_exp5, count, ancilla_exp5 = recursive_karatsuba(eng, t1_exp5, t2_copy, n, count, ancilla_exp5)
    Reduction(eng, t1_exp5)

    t1_exp5 = Squaring_temp(eng, t1_exp5)
    count = 0
    out_exp5, count, ancilla_exp5 = recursive_karatsuba(eng, t1_exp5, input, n, count, ancilla_exp5)
    Reduction(eng, out_exp5)

    ###### exp_7 is finished ######

    ###### exp_29 start ###########

    count = 0
    t2_exp29, count, ancilla_exp29, = recursive_karatsuba(eng, t1_exp29, input, n, count, ancilla_exp29)
    Reduction(eng, t2_exp29)

    t2_exp29_copy = eng.allocate_qureg(n)
    copy(eng, t2_exp29, t2_exp29_copy, n)

    t1_exp29 = Squaring_temp(eng, t2_exp29)
    t1_exp29 = Squaring_temp(eng, t1_exp29)
    t1_exp29 = Squaring_temp(eng, t1_exp29)

    count = 0
    t2_exp29, count, ancilla_exp29, = recursive_karatsuba(eng, t1_exp29, t2_exp29_copy, n, count, ancilla_exp29)
    Reduction(eng, t2_exp29)

    t1_exp29 = Squaring_temp(eng, t2_exp29)
    count = 0
    t2_exp29, count, ancilla_exp29, = recursive_karatsuba(eng, t1_exp29, input, n, count, ancilla_exp29)
    Reduction(eng, t2_exp29)

    t2_exp29_copy = eng.allocate_qureg(n)
    copy(eng, t2_exp29, t2_exp29_copy, n)

    t1_exp29 = Squaring_temp(eng, t2_exp29)
    for i in range(6):
        t1_exp29 = Squaring_temp(eng, t1_exp29)

    count = 0
    t2_exp29, count, ancilla_exp29, = recursive_karatsuba(eng, t1_exp29, t2_exp29_copy, n, count, ancilla_exp29)
    Reduction(eng, t2_exp29)

    ##
    t2_exp29_copy = eng.allocate_qureg(n)
    copy(eng, t2_exp29, t2_exp29_copy, n)

    t1_exp29 = Squaring_temp(eng, t2_exp29)
    for i in range(13):
        t1_exp29 = Squaring_temp(eng, t1_exp29)

    count = 0
    t2_exp29, count, ancilla_exp29, = recursive_karatsuba(eng, t1_exp29, t2_exp29_copy, n, count, ancilla_exp29)
    Reduction(eng, t2_exp29)

    t1_exp29 = Squaring_temp(eng, t2_exp29)
    count = 0
    out_exp29, count, ancilla_exp29, = recursive_karatsuba(eng, t1_exp29, input, n, count, ancilla_exp29)
    Reduction(eng, out_exp29)

    return out_exp5, out_exp29, ancilla_exp5

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
    temp = eng.allocate_qureg(3)

    for i in range(96):
        new_in.append(input[i])
        # print(count,'(',count2,')')
        new_in.append(input[96 + i])
        # print(count + 96, '(', count2, ')')

    CNOT | (input[189], temp[0])
    CNOT | (input[190], temp[1])
    CNOT | (input[191], temp[2])

    for i in reversed(range(93)):
        CNOT | (input[96 + i], new_in[7 + 2 * i])
        CNOT | (input[96 + i], new_in[2 + 2 * i])
        CNOT | (input[96 + i], new_in[0 + 2 * i])

    CNOT | (temp[0], new_in[8])
    CNOT | (temp[0], new_in[3])
    CNOT | (temp[0], new_in[2])
    CNOT | (temp[0], new_in[1])

    CNOT | (temp[0], new_in[188])
    CNOT | (temp[0], new_in[186])

    CNOT | (temp[1], new_in[10])
    CNOT | (temp[1], new_in[5])
    CNOT | (temp[1], new_in[4])
    CNOT | (temp[1], new_in[3])

    CNOT | (temp[1], new_in[190])
    CNOT | (temp[1], new_in[188])

    CNOT | (temp[2], new_in[190])

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

def Reduction(eng, x):

    for i in range(10):
        CNOT | (x[192 + i], x[0 + i])
        CNOT | (x[192 + i], x[1 + i])
        CNOT | (x[192 + i], x[2 + i])
        CNOT | (x[192 + i], x[7 + i])

        CNOT | (x[202 + i], x[10 + i])
        CNOT | (x[202 + i], x[11 + i])
        CNOT | (x[202 + i], x[12 + i])
        CNOT | (x[202 + i], x[17 + i])

        CNOT | (x[212 + i], x[20 + i])
        CNOT | (x[212 + i], x[21 + i])
        CNOT | (x[212 + i], x[22 + i])
        CNOT | (x[212 + i], x[27 + i])

        CNOT | (x[222 + i], x[30 + i])
        CNOT | (x[222 + i], x[31 + i])
        CNOT | (x[222 + i], x[32 + i])
        CNOT | (x[222 + i], x[37 + i])

        CNOT | (x[232 + i], x[40 + i])
        CNOT | (x[232 + i], x[41 + i])
        CNOT | (x[232 + i], x[42 + i])
        CNOT | (x[232 + i], x[47 + i])

        CNOT | (x[242 + i], x[50 + i])
        CNOT | (x[242 + i], x[51 + i])
        CNOT | (x[242 + i], x[52 + i])
        CNOT | (x[242 + i], x[57 + i])

        CNOT | (x[252 + i], x[60 + i])
        CNOT | (x[252 + i], x[61 + i])
        CNOT | (x[252 + i], x[62 + i])
        CNOT | (x[252 + i], x[67 + i])

        CNOT | (x[262 + i], x[70 + i])
        CNOT | (x[262 + i], x[71 + i])
        CNOT | (x[262 + i], x[72 + i])
        CNOT | (x[262 + i], x[77 + i])

        CNOT | (x[272 + i], x[80 + i])
        CNOT | (x[272 + i], x[81 + i])
        CNOT | (x[272 + i], x[82 + i])
        CNOT | (x[272 + i], x[87 + i])

        CNOT | (x[282 + i], x[90 + i])
        CNOT | (x[282 + i], x[91 + i])
        CNOT | (x[282 + i], x[92 + i])
        CNOT | (x[282 + i], x[97 + i])

        CNOT | (x[292 + i], x[100 + i])
        CNOT | (x[292 + i], x[101 + i])
        CNOT | (x[292 + i], x[102 + i])
        CNOT | (x[292 + i], x[107 + i])

        CNOT | (x[302 + i], x[110 + i])
        CNOT | (x[302 + i], x[111 + i])
        CNOT | (x[302 + i], x[112 + i])
        CNOT | (x[302 + i], x[117 + i])

        CNOT | (x[312 + i], x[120 + i])
        CNOT | (x[312 + i], x[121 + i])
        CNOT | (x[312 + i], x[122 + i])
        CNOT | (x[312 + i], x[127 + i])

        CNOT | (x[322 + i], x[130 + i])
        CNOT | (x[322 + i], x[131 + i])
        CNOT | (x[322 + i], x[132 + i])
        CNOT | (x[322 + i], x[137 + i])

        CNOT | (x[332 + i], x[140 + i])
        CNOT | (x[332 + i], x[141 + i])
        CNOT | (x[332 + i], x[142 + i])
        CNOT | (x[332 + i], x[147 + i])

        CNOT | (x[342 + i], x[150 + i])
        CNOT | (x[342 + i], x[151 + i])
        CNOT | (x[342 + i], x[152 + i])
        CNOT | (x[342 + i], x[157 + i])

        CNOT | (x[352 + i], x[160 + i])
        CNOT | (x[352 + i], x[161 + i])
        CNOT | (x[352 + i], x[162 + i])
        CNOT | (x[352 + i], x[167 + i])

        CNOT | (x[362 + i], x[170 + i])
        CNOT | (x[362 + i], x[171 + i])
        CNOT | (x[362 + i], x[172 + i])
        CNOT | (x[362 + i], x[177 + i])

    for i in range(5):
        CNOT | (x[372 + i], x[180 + i])
        CNOT | (x[372 + i], x[181 + i])
        CNOT | (x[372 + i], x[182 + i])
        CNOT | (x[372 + i], x[187 + i])

    #[249]
    CNOT | (x[377], x[185])
    CNOT | (x[377], x[186])
    CNOT | (x[377], x[187])

    CNOT | (x[378], x[1])
    CNOT | (x[378], x[2])
    CNOT | (x[378], x[3])
    CNOT | (x[378], x[8])

    CNOT | (x[377], x[0])
    CNOT | (x[377], x[1])
    CNOT | (x[377], x[2])
    CNOT | (x[377], x[7])

    # [250]
    CNOT | (x[378], x[186])
    CNOT | (x[378], x[187])
    CNOT | (x[378], x[188])

    # [251]
    CNOT | (x[379], x[187])
    CNOT | (x[379], x[188])
    CNOT | (x[379], x[189])

    CNOT | (x[380], x[3])
    CNOT | (x[380], x[4])
    CNOT | (x[380], x[5])
    CNOT | (x[380], x[10])

    CNOT | (x[379], x[2])
    CNOT | (x[379], x[3])
    CNOT | (x[379], x[4])
    CNOT | (x[379], x[9])

    # [252]
    CNOT | (x[380], x[188])
    CNOT | (x[380], x[189])
    CNOT | (x[380], x[190])

    # [253]
    CNOT | (x[381], x[189])
    CNOT | (x[381], x[190])
    CNOT | (x[381], x[191])

    CNOT | (x[382], x[0])
    CNOT | (x[382], x[1])
    CNOT | (x[382], x[2])

    CNOT | (x[382], x[5])
    CNOT | (x[382], x[6])
    CNOT | (x[382], x[12])

    CNOT | (x[381], x[4])
    CNOT | (x[381], x[5])
    CNOT | (x[381], x[6])
    CNOT | (x[381], x[11])

    # [253]
    CNOT | (x[382], x[190])
    CNOT | (x[382], x[191])

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

print('\n')
print('Estimate cost...')
Resource = ResourceCounter()
eng = MainEngine(Resource)
resource_check = 1
AIM(eng)

print(Resource)
eng.flush()
