# Elliot Bidwell
# CS 2300 - 001
# Programming Project 2

from math import *


class LinearSystem:

    def __init__(self, matrix):
        self.matrix = matrix

    def set_matrix(self, matrix):
        self.matrix = matrix

    # This function employs Cramer's Rule to calculate the solutions of a linear system
    def cramer(self, matrix, vector):
        vect_a1 = matrix[0]
        vect_a2 = matrix[1]
        vect_b = vector

        determinant_of_A = self.determinant_2D(matrix)

        x1 = (vect_b[0] * vect_a2[1] - vect_b[1] * vect_a2[0]) / determinant_of_A
        x2 = (vect_a1[0] * vect_b[1] - vect_a1[1] * vect_b[0]) / determinant_of_A

        return [x1, x2]

    # This function calculates the determinant of a 2x2 matrix
    def determinant_2D(self, matrix):
        return matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]

    # This function returns the result of multiplying two 2D arrays. It multiplies its second argument by its second
    def matrix_mult2D(self, matrix1, matrix2):
        return [[matrix1[0][0] * matrix2[0][0] + matrix1[1][0] * matrix2[0][1],
                 matrix1[0][1] * matrix2[0][0] + matrix1[1][1] * matrix2[0][1]],
                [matrix1[0][0] * matrix2[1][0] + matrix1[1][0] * matrix2[1][1],
                 matrix1[0][1] * matrix2[1][0] + matrix1[1][1] * matrix2[1][1]]]

    # This function solves a 2x2 linear system and returns the two solutions in a 2x1 matrix or list.
    # It checks to see if the determinant of A = 0, or that A is linearly dependant, and uses Cramer's Rule
    # if it is not. Otherwise, it decides whether the system is undetermined or unsolvable.
    def solveSystem(self):
        matrix_A = [self.matrix[0], self.matrix[1]]
        vect_b = self.matrix[2]

        if self.determinant_2D(matrix_A) == 0:
            if self.determinant_2D([matrix_A[0], vect_b]) == 0:
                result = f'System underdetermined'
            else:
                result = f'System inconsistent'
        else:
            result = self.cramer(matrix_A, vect_b)
        return result

    def eigenThings(self):
        matrix_A = [self.matrix[0], self.matrix[1]]

        # Calculate the coefficients of the quadratic equation
        # First coefficient is always 1
        a = 1
        # Calculate 2nd coefficient
        b = (-matrix_A[0][0]) + (-matrix_A[1][1])
        # Calculate 3rd coefficient
        c = (matrix_A[0][0] * matrix_A[1][1]) - (matrix_A[0][1] * matrix_A[1][0])

        # Exception handler for the case that calculating either eigenvalue results in no real solution (i.e. the square
        # root of a negative number)
        try:
            # Use the quadratic equation to find the eigenvalues
            lambda1 = ((-b) + sqrt((b ** 2) - 4 * a * c))/ 2 * a
            lambda2 = ((-b) - sqrt((b ** 2) - 4 * a * c)) / 2 * a

            # Check which lambda has the higher absolute value and swap them accordingly
            if abs(lambda2) > abs(lambda1):
                lambda1, lambda2 = lambda2, lambda1

            # Compose the matrix for the homogeneous system for the dominant eigenvector
            r1_matrix = [[matrix_A[0][0] - lambda1, matrix_A[0][1]], [matrix_A[1][0], matrix_A[1][1] - lambda1]]

            # Compose the matrix for the homogeneous system for the non-dominant eigenvector
            r2_matrix = [[matrix_A[0][0] - lambda2, matrix_A[0][1]], [matrix_A[1][0], matrix_A[1][1] - lambda2]]

            # Create row and column pivoted matrices for r1 and r2. NOT USED
            row_pivoted_r1 = [[r1_matrix[0][1], r1_matrix[0][0]], [r1_matrix[1][1], r1_matrix[1][0]]]
            row_pivoted_r2 = [[r2_matrix[0][1], r2_matrix[0][0]], [r2_matrix[1][1], r2_matrix[1][0]]]
            column_pivoted_r1 = [r1_matrix[1], r1_matrix[0]]
            column_pivoted_r2 = [r2_matrix[1], r2_matrix[0]]

            if r1_matrix[0][1] != 0:

                # Compose a shear matrix to do forward elimination with the r1_matrix
                shear_matrix_r1 = [[1, -(r1_matrix[0][1] / r1_matrix[0][0])], [0, 1]]
                # Shear r1_matrix via matrix multiplication
                sheared_r1 = self.matrix_mult2D(shear_matrix_r1, r1_matrix)

                # Set r12 to 1 and use back substitution to solve for r1
                r12 = 1
                r11 = -(sheared_r1[1][0]) / sheared_r1[0][0]
                r1_magnitude = sqrt(r11 ** 2 + r12 ** 2)
                r11 = r11 / r1_magnitude
                r12 = r12 / r1_magnitude
                r1 = [r11, r12]

            if r2_matrix[0][1] != 0:

                # Compose a shear matrix to do forward elimination with the r2_matrix
                shear_matrix_r2 = [[1, -(r2_matrix[0][1] / r2_matrix[0][0])], [0, 1]]
                # Shear r2_matrix via matrix multiplication
                sheared_r2 = self.matrix_mult2D(shear_matrix_r2, r2_matrix)

                # Set r22 to 1 and use back substitution to solve for r2
                r22 = 1
                r21 = -(sheared_r2[1][0]) / sheared_r2[0][0]
                r2_magnitude = sqrt(r21 ** 2 + r22 ** 2)
                r21 = r21 / r2_magnitude
                r22 = r22 / r2_magnitude
                r2 = [r21, r22]

            # Compose a matrix containing the eigenvectors
            matrix_R = [r1, r2]

            # Compose a matrix containing the eigenvalues
            lambda_matrix = [[lambda1, float(0)], [float(0), lambda2]]

            # Compose the transpose of R
            R_transpose_tuples = list(zip(*matrix_R))
            matrix_R_transpose = [list(column) for column in R_transpose_tuples]

            result = [[matrix_R, lambda_matrix, matrix_R_transpose], int((self.matrix_mult2D(self.matrix_mult2D(matrix_R, lambda_matrix), matrix_R_transpose)) == matrix_A)]

        except ValueError:
            result = f'No real eigenvalues.'

        return result

    def area_distance(self):
        p1 = self.matrix[0]
        p2 = self.matrix[1]
        p3 = self.matrix[2]
        v1 = [(p2[i] - p1[i]) for i in range(len(p1))]
        v2 = [(p3[i] - p1[i]) for i in range(len(p1))]

        if len(v1) == 2:
            area = self.determinant_2D([v1, v2]) / 2

            # Create a
            a, b = -(v1[1]), v1[0]
            c = -(a * p1[0])-(b * p1[1])
            distance = a * p3[0] + b * p3[1] + c

        else:
            # Compute the cross product of the vectors that span the parallelogram
            cross_product = [self.determinant_2D([[v1[1], v1[2]], [v2[1], v2[2]]]),
                             -(self.determinant_2D([[v1[0], v1[2]], [v2[0], v2[2]]])),
                             self.determinant_2D([[v1[0], v1[1]], [v2[0], v2[1]]])]
            # Compute the magnitude of the cross product and divide it by 2 to get the area of the triangle
            area = sqrt(cross_product[0] ** 2 + cross_product[1] ** 2 + cross_product[2] ** 2) / 2

            # Compute the magnitude of the vector that passes between the two bisected points
            v1_magnitude = sqrt(v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2)

            # Populate a normalized vector of v1
            normal = [(component/v1_magnitude) for component in v1]

            # Set the first 3 coefficients of the plane equation to the components of the normal vector
            a, b, c = normal
            # Calculate a point on the plane by finding the midpoint of the tow bisected points
            p = [(p1[0] + p2[0])/2, (p1[1] + p2[1])/2, (p1[2] + p2[2])/2]

            # Calculate the fourth coefficient using the midpoint
            d = -(a * p[0] + b * p[1] + c * p[2])

            # Calculate the distance by plugging p3 into the plane equation
            x1, x2, x3 = p3
            distance = ((a * x1) + (b * x2) + (c * x3) + d)
        return round(area, 4), round(distance, 4)


# This function prints the solution to the 2x3 system of equations
def printSysSolution(solution, filename):
    output_file = open(filename, 'w', encoding="utf-8")
    if isinstance(solution, str):
        output_file.write(solution)
    else:
        output_file.write(f'{solution[0]}\n{solution[1]}')
    output_file.close()


# This functions prints the result of calculating eigen things
def printEigen(eigen_things, filename):
    output_file = open(filename, 'w', encoding="utf-8")

    if isinstance(eigen_things, str):
        output_file.write(eigen_things)

    else:
        is_decomposition_equal = eigen_things[1]
        eigen_things = eigen_things[0]
        transpose = list(zip(*eigen_things[1]))
        lambda_matrix = [list(column) for column in transpose]
        for row in lambda_matrix:
            for element in row:
                output_file.write(f'{round(element, 4)}\t')
            output_file.write('\n')

        transpose = list(zip(*eigen_things[0]))
        matrix_R = [list(column) for column in transpose]
        for row in matrix_R:
            for element in row:
                output_file.write(f'{round(element, 4)}\t')
            output_file.write('\n')

        transpose = list(zip(*eigen_things[2]))
        matrix_R_transpose = [list(column) for column in transpose]

        eigendecomposition = [matrix_R, lambda_matrix, matrix_R_transpose]
        for matrix in eigendecomposition:
            for element in matrix[0]:
                output_file.write(f'{round(element, 4)}\t')
        output_file.write('\n')

        for matrix in eigendecomposition:
            for element in matrix[1]:
                output_file.write(f'{round(element, 4)}\t')
        output_file.write('\n')

        output_file.write(f'{is_decomposition_equal}')

    output_file.close()


# This function prints the results of calculating the area of a triangle and the distance from a plane/line
def print_area_distance(area_distance, filename):
    output_file = open(filename, 'w', encoding="utf-8")
    area, distance = area_distance
    output_file.write(f'{area}\n{distance}')
    output_file.close()


# This file receives input from the input file
def get_input(filename):
    # Create a new file object
    input_file = open(filename, 'r', encoding="utf-8")
    input_values = [[int(x) for x in line.split()] for line in input_file]
    input_file.close()

    # The file is parsed by row and stores the elements of each row in one of the lists contained within
    # input_values. This call to zip() makes a list of tuples that is a transpose of input_values and
    # contains the columns of the matrix in the file.
    transpose = list(zip(*input_values))

    # Cast each tuple or column in transpose to a list and populate input_values with them
    input_values = [list(column) for column in transpose]

    return input_values


input_matrix = get_input('2x3_input1.txt')
system = LinearSystem(input_matrix)
printSysSolution(system.solveSystem(), '2x3_system_output1.txt')
printEigen(system.eigenThings(), '2x3_eigen_output1.txt')
print_area_distance(system.area_distance(), '2x3_output1.txt')

input_matrix = get_input('2x3_input2.txt')
system.set_matrix(input_matrix)
printSysSolution(system.solveSystem(), '2x3_system_output2.txt')
printEigen(system.eigenThings(), '2x3_eigen_output2.txt')
print_area_distance(system.area_distance(), '2x3_output2.txt')

input_matrix = get_input('2x3_input_rotation.txt')
system.set_matrix(input_matrix)
printSysSolution(system.solveSystem(), "2x3_system_output_rotation.txt")
printEigen(system.eigenThings(), '2x3_eigen_output_rotation.txt')
print_area_distance(system.area_distance(), '2x3_area_distance_output_rotation.txt')

input_matrix = get_input('2x3_input_pivot.txt')
system.set_matrix(input_matrix)
printSysSolution(system.solveSystem(), '2x3_system_output_pivot')
printEigen(system.eigenThings(), '2x3_eigen_output_pivot')
print_area_distance(system.area_distance(), '2x3_area_distance_output_pivot.txt')

input_matrix = get_input('2x3_input_lin_dependent.txt')
system.set_matrix(input_matrix)
printSysSolution(system.solveSystem(), '2x3_system_output_lin_dependent.txt')
printEigen(system.eigenThings(), '2x3_eigen_output_lin_dependent.txt')
print_area_distance(system.area_distance(), '2x3_area_distance_output_lin_dependent.txt')

input_matrix = get_input('3x3_input1.txt')
system.set_matrix(input_matrix)
print_area_distance(system.area_distance(), '3x3_output1.txt')

input_matrix = get_input('3x3_input2.txt')
system.set_matrix(input_matrix)
print_area_distance(system.area_distance(), '3x3_output2.txt')


