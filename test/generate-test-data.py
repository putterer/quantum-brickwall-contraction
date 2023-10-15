import numpy as np

iterations = 1

qubits = 8
statevec_size = 2**qubits
gate = np.random.uniform(-1.0, 1.0, (4, 4)) + 1j * np.random.uniform(-1.0, 1.0, (4, 4))

# one-hot vector
input = np.random.uniform(-1.0, 1.0, (statevec_size)) + 1j * np.random.uniform(-1.0, 1.0, (statevec_size))

# operation
op = np.ones((1, 1), dtype=complex)
for i in range(0, qubits, 2):
    op = np.kron(op, gate)

# apply
output = input
for i in range(iterations):
    output = op @ output


data = f"TestData<{qubits}> TEST_GENERATED = {'{'}\n"

# parameters
data += f"    .iterations = {iterations},\n"

# input
data += "    .input = {"
for entry in input:
    data += f"complex({entry.real},{entry.imag})"
    data += ","
data = data[0:-1]
data += "},\n"

# output
data += "    .output = {"
for entry in output:
    data += f"complex({entry.real},{entry.imag})"
    data += ","
data = data[0:-1]
data += "},\n"

# gate
data += "    .gate = {"
for row in range(gate.shape[0]):
    data += "{"
    for column in range(gate.shape[1]):
        entry = gate[row][column]
        data += f"complex({entry.real},{entry.imag})"
        data += ","
    data = data[0:-1]
    data += "},"
data = data[0:-1]
data += "}\n"



data += "};"

print(data)

