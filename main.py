import numpy as np

def to_row_echelon_form(matrix):
    matrix = matrix.astype(float)
    rows, cols = matrix.shape
    lead = 0

    for r in range(rows):
        if lead >= cols:
            return matrix

        i = r
        while matrix[i, lead] == 0:
            i += 1
            if i == rows:
                i = r
                lead += 1
                if cols == lead:
                    return matrix

        matrix[[i, r]] = matrix[[r, i]]

        lv = matrix[r, lead]
        matrix[r] = matrix[r] / lv

        for i in range(r + 1, rows):
            matrix[i] = matrix[i] - matrix[r] * matrix[i, lead]

        lead += 1

    return matrix

def matrix_to_equations(matrix):
    rows, cols = matrix.shape
    equations = []

    for r in range(rows):
        equation = []
        if np.all(matrix[r, :-1] == 0):
            continue

        for c in range(cols - 1):
            coeff = matrix[r, c]
            variable = f"x_{c + 1}"

            if coeff > 0:
                if equation:
                    equation.append(f"+ {coeff}*{variable}" if coeff != 1 else f"+ {variable}")
                else:
                    equation.append(f"{coeff}*{variable}" if coeff != 1 else f"{variable}")
            elif coeff < 0:
                equation.append(f"- {-coeff}*{variable}" if coeff != -1 else f"- {variable}")

        free_term = matrix[r, -1]
        equation_str = " ".join(equation) + f" = {free_term}"
        equations.append(equation_str)

    return equations

def calculate_rank(matrix):
    rank = 0
    for row in matrix:
        if np.any(row != 0):
            rank += 1
    return rank

def calculate_one_solution(matrix):
    rows, cols = matrix.shape
    free_terms = matrix[:, -1]
    variables = np.zeros(rows)

    for i in range(rows - 1, -1, -1):
        variables[i] = free_terms[i]

        for j in range(i + 1, cols - 1):
            variables[i] -= matrix[i, j] * variables[j]

        if matrix[i, i] != 0:
            variables[i] /= matrix[i, i]

    return variables

def calculate_many_solutions(matrix):
    num_rows, num_cols = matrix.shape

    expressions = ['0'] * (num_cols - 1)

    for i in range(num_rows - 1, -1, -1):
        if np.any(matrix[i, :-1]):
            for j in range(num_cols - 1):
                if matrix[i, j] == 1:
                    right_side = matrix[i, -1]
                    expr = str(right_side)

                    for k in range(num_cols - 1):
                        if matrix[i, k] != 0 and k != j:
                            coeff = matrix[i, k]
                            if coeff > 0:
                                expr += f' - {coeff}*x{k + 1}'
                            else:
                                expr += f' + {abs(coeff)}*x{k + 1}'

                    expressions[j] = expr.strip()

    for i, expr in enumerate(expressions):
        if expr != '0':
            print(f'x{i + 1} = {expr}')

def calculate_answer(matrix, matrix_aug):
    rank_matrix = calculate_rank(to_row_echelon_form(matrix))
    rank_matrix_aug = calculate_rank(to_row_echelon_form(matrix_aug))
    n = len(matrix[0])

    if rank_matrix < rank_matrix_aug: # НЕТ РЕШЕНИЙ
        print("Система не имеет решений!")
        return

    if rank_matrix < n: # БЕСКОНЕЧНО РЕШЕНИЙ

        print("Множество решений")

        calculate_many_solutions(to_row_echelon_form(matrix_aug))
        return

    # ОДНО РЕШЕНИЕ
    solution = calculate_one_solution(to_row_echelon_form(matrix_aug))

    print("Единственное решение:")
    for i, value in enumerate(solution, start=1):
        print(f"x_{i} = {value:.2f}")
    return

if __name__ == "__main__":

    matrix_aug4 = np.array([[7, -2, -1, 2],
                       [6, -4, -5, 3],
                       [1, 2, 4, 5]])
    matrix4 = matrix_aug4[:, :-1]

    matrix_aug2 = np.array([[4, 3, 3, 2],
                       [-5, 4, 3, 2],
                       [2, 3, 4, 8]])
    matrix2 = matrix_aug2[:, :-1]

    matrix_aug3 = np.array([[2, 3, -1, 1, 1],
                       [8, 12, -9, 8, 3],
                       [4, 6, 3, -2, 3],
                       [2, 3, 9, -7, 3]])
    matrix3 = matrix_aug3[:, :-1]

    # Ввод матрицы (расширенной)
    matrix_aug = np.array([[2, -1, 3, -5, 1],
                       [1, -1, -5, 0, 2],
                       [3, -2, -2, -5, 3],
                       [7, -5, -9, -10, 8]])
    print("Расширенная матрица:")
    print(matrix_aug)

    # Обычная матрица
    matrix = matrix_aug[:, :-1]
    print("Изначальная матрица:")
    print(matrix)

    # Приведение к ступенчатому виду
    row_echelon_matrix = to_row_echelon_form(matrix_aug4)
    print("Ступенчатый вид матрицы:")
    print(row_echelon_matrix)

    calculate_answer(matrix4, matrix_aug4)
