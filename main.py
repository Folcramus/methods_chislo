import numpy as np


def to_row_echelon_form(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Приведение расширенной матрицы системы к ступенчатому виду."""
    matrix = np.hstack((a, b.reshape(-1, 1))).astype(float)  # Приводим к float для точности
    rows, cols = matrix.shape
    lead = 0  # Ведущий элемент
    print("Матрица в ступенчатом виде:")
    for r in range(rows):
        if lead >= cols:
            print(matrix)
            return matrix

        i = r
        while np.isclose(matrix[i, lead], 0):
            i += 1
            if i == rows:
                i = r
                lead += 1
                if lead == cols:
                    print(matrix)
                    return matrix

        # Меняем местами строки
        matrix[[i, r]] = matrix[[r, i]]

        # Нормализуем ведущий элемент
        matrix[r] = matrix[r] / matrix[r, lead]

        # Обнуляем элементы ниже ведущего
        for i in range(r + 1, rows):
            matrix[i] -= matrix[r] * matrix[i, lead]

        lead += 1
    print(matrix)
    return matrix


def back_substitution(matrix: np.ndarray) -> np.ndarray:
    """Обратный ход Гаусса для приведения матрицы к верхнетреугольному виду."""
    rows, cols = matrix.shape
    cols -= 1  # Последний столбец — это столбец свободных коэффициентов

    for i in range(rows - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            factor = matrix[j, i]
            matrix[j] -= factor * matrix[i]

    return matrix


def solving_equations(matrix: np.ndarray) -> None:
    """Решение системы уравнений на основе приведенной к ступенчатому виду матрицы."""
    rank_a = np.linalg.matrix_rank(matrix[:, :-1])
    rank_aug = np.linalg.matrix_rank(matrix)

    if rank_a < rank_aug:
        print("Система не имеет решений")
        print('-------------------------------------------------')
        return

    # Обратный ход (нахождение решения)

    back_substitution(matrix)
    print("Решенная система уравнений:")
    print(matrix)
    n = len(matrix[0])
    print()

    # Параметрическое решение
    if rank_a < n:
        calculate_many_solutions(matrix)
        print('-------------------------------------------------')
    else:
        solution = calculate_one_solution(matrix)
        print("Единственное решение:")
        for i, value in enumerate(solution, start=1):
            print(f"x_{i} = {value:.8f}")
        print('-------------------------------------------------')


def calculate_one_solution(matrix: np.ndarray) -> np.ndarray:
    """Вычисление единственного решения СЛАУ методом обратного хода."""
    rows, cols = matrix.shape
    free_terms = matrix[:, -1]
    variables = np.zeros(rows)

    for i in range(rows - 1, -1, -1):
        variables[i] = free_terms[i] - np.dot(matrix[i, i + 1:cols - 1], variables[i + 1:])
        if matrix[i, i] != 0:
            variables[i] /= matrix[i, i]

    return variables


def calculate_many_solutions(matrix: np.ndarray) -> None:
    """Вывод параметрических решений для недоопределенных систем."""
    num_rows, num_cols = matrix.shape
    expressions = ['0'] * (num_cols - 1)

    for i in range(num_rows - 1, -1, -1):
        if np.any(matrix[i, :-1]):
            for j in range(num_cols - 1):
                if matrix[i, j] == 1:
                    expr = f"{matrix[i, -1]:.8f}"
                    for k in range(num_cols - 1):
                        coeff = matrix[i, k]
                        if coeff != 0 and k != j:
                            sign = '-' if coeff > 0 else '+'
                            expr += f" {sign} {abs(coeff)}*x{k + 1}"
                    expressions[j] = expr.strip()

    for i, expr in enumerate(expressions):
        if expr != '0':
            print(f"x{i + 1} = {expr}")


if __name__ == "__main__":
    A = np.array([[2, -1, 3, -5],
                  [1, -1, -5, 0],
                  [3, -2, -2, -5],
                  [7, -5, -9, -10]], dtype=float)

    b = np.array([1, 2, 3, 8], dtype=float)

    A1 = np.array([[7, -2, -1],
                   [6, -4, -5],
                   [1, 2, 4]], dtype=float)



    b1 = np.array([2, 3, 5], dtype=float)

    A2 = np.array([[9, 8, 7],
                   [6, 5, 4],
                   [3, 2, 1]])
    b2 = np.array([8, 5, 3])

    # меняй тут
    Ag = np.array([[0.87, -0.27, 0.22, 0.18],
                   [-0.21, 1, 0.45, -0.18],
                   [-0.12, -0.13, 1.33, -0.18],
                   [-0.33, 0.05, -0.06, 1.28]])
    bg = np.array([1.21, -0.33, -0.48, -0.17])
    solving_equations(to_row_echelon_form(A, b))
    solving_equations(to_row_echelon_form(A1, b1))
    solving_equations(to_row_echelon_form(A2, b2))
    solving_equations(to_row_echelon_form(Ag, bg))
