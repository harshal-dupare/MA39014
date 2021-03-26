import numpy as np

"""

18MA20015 | Harshal Dupare
Lab 7 | 24-2-2021

"""

class Revised_Simplex:
    def __init__(self):
        self.WIDTH = 9
        self.MAX_itr = 25
        self.n = None
        self.m = None
        self.C_j = None
        self.A_ij = None
        self.X_B = None
        self.B = None
        self.C_b = None
        self.B_inv = None
        self.y = None
        self.basic_var = None
        self.non_basic_var = None
        pass

    def file_input(self, file_name):
        question_file = open(file_name)

        self.n = int(question_file.readline())
        self.m = int(question_file.readline())

        ls = (question_file.readline()).split()

        # maximize the function
        self.C_j = np.zeros(self.n+self.m)
        for i in range(self.n):
            self.C_j[i] = float(ls[i])

        self.A_ij = np.zeros((self.m, self.m+self.n))
        self.X_B = np.zeros(self.m)
        # all "<=" type equations
        for i in range(self.m):
            ls = (question_file.readline()).split()
            for j in range(self.n):
                self.A_ij[i, j] = float(ls[j])
            self.A_ij[i, (self.n+i)] = 1
            self.X_B[i] = ls[self.n]

        question_file.close()
        pass
    
    def input(self):
        print("All conditions must be of type <= and problem must be of type Maximization\n")
        self.n = int(input("Enter number of variables: "))
        self.m = int(input("Enter number of equations: "))

        print('Enter the coefficient of x_i\'s for Z in format\nC_1 C_2 ... C_n')
        ls = (input()).split()

        # maximize the function
        self.C_j = np.zeros(self.n+self.m)
        for i in range(self.n):
            self.C_j[i] = float(ls[i])

        self.A_ij = np.zeros((self.m, self.m+self.n))
        self.X_B = np.zeros(self.m)

        # all "<=" type equations
        print(f"For next {self.m} lines \nEnter the coefficient in the format \"a_i1 a_i2  ... a_in b_i\" in equation ith for the ith line\n")
        for i in range(self.m):
            print(f"Enter the {i+1}'th equation: ")
            ls = (input()).split()
            for j in range(self.n):
                self.A_ij[i, j] = float(ls[j])
            self.A_ij[i, (self.n+i)] = 1
            self.X_B[i] = ls[self.n]

        pass

    def printtable(self):
        print(" ||   Initial Table is  ||")
        print(f"{'C_b':<{self.WIDTH}}",end="")
        print(f"{'X_b':<{self.WIDTH}}",end="")
        for i in range(self.n+self.m):
            print(f"{'x_'+str(i+1):<{self.WIDTH}}",end="")
        print(f"{'b':<{self.WIDTH}}")

        for i in range(self.m):
            print(f"{self.C_b[i]:<{self.WIDTH}}",end="")
            print(f"{'x_'+str(self.basic_var[i]+1):<{self.WIDTH}}",end="")
            for j in range(self.n+self.m):
                print(f"{self.A_ij[i][j]:<{self.WIDTH}}",end="")
            print(f"{self.X_B[i]:<{self.WIDTH}}")

        pass

    def initilize(self):

        print('\n\nX_B : ' + str(self.X_B))

        self.B = np.eye(self.m)  
        print('\nB : \n'+str(self.B))

        self.C_b = np.zeros(self.m)  
        print('\nC_b for the basic varaibles: ' + str(self.C_b))
        
        self.B_inv = np.linalg.inv(self.B)  
        print('\nB_inverse : \n' + str(self.B_inv))

        self.y = np.dot(self.C_b, self.B_inv)  
        print('\ny = C_b*B_inverse = ' + str(self.y))

        self.basic_var = np.arange(self.m)+self.n
        self.non_basic_var = np.arange(self.n)

        bsvrs = ['x_'+str(bvi+1) for bvi in self.basic_var]
        print('\nBasic variables: ' + str(bsvrs))
        nbsvrs = ['x_'+str(bvi+1) for bvi in self.non_basic_var]
        print('Non-Basic variables: ' + str(nbsvrs))

        start_id = 0
        non_basic_coeff_id = 0
        max_val = 0
        for i in range(self.n):
            if (self.C_j[self.non_basic_var[i]] - np.dot(self.y, np.transpose(self.A_ij[:, self.non_basic_var[i]]))) > max_val:
                max_val = self.C_j[self.non_basic_var[i]] - np.dot(self.y, np.transpose(self.A_ij[:, self.non_basic_var[i]]))
                start_id = self.non_basic_var[i]
                non_basic_coeff_id = i

        print(f'\n\nMax_j Cj - Zj = {max_val}')
        print(f'\nEntering variable: x_{start_id+1}')

        return start_id, non_basic_coeff_id

    def compute(self, start_id, non_basic_coeff_id):

        iter_number = 0
        next_loop = True

        while next_loop :
            print(f"\n\n||   Iteration number : {iter_number+1}   ||")
            next_loop = False

            arr = np.divide(self.X_B, np.dot(self.B_inv,  self.A_ij[:, start_id]))

            iter_id = np.where(arr > 0, arr, np.inf).argmin()
            minim = np.where(arr > 0, arr, np.inf).min()

            self.non_basic_var[non_basic_coeff_id] = self.basic_var[iter_id]
            self.basic_var[iter_id] = start_id

            self.X_B = (self.X_B) - minim * np.dot(self.B_inv,  self.A_ij[:, start_id])
            self.X_B[iter_id] = minim
            print('\n\nX_B : ' + str(self.X_B))

            self.C_b[iter_id] = self.C_j[self.basic_var[iter_id]]
            print('\nC_b for the basic varaibles: ' + str(self.C_b))

            self.B[:, iter_id] = self.A_ij[:, start_id]
            print('\nB : \n'+str(self.B))

            self.B_inv = np.linalg.inv(self.B)
            print('\nB_inverse : \n' + str(self.B_inv))

            self.y = np.dot(self.C_b, self.B_inv)
            print('\ny = C_b*B_inverse = ' + str(self.y))

            bsvrs = ['x_'+str(bvi+1) for bvi in self.basic_var]
            print('\nBasic variables: ' + str(bsvrs))
            nbsvrs = ['x_'+str(bvi+1) for bvi in self.non_basic_var]
            print('Non-Basic variables: ' + str(nbsvrs))

            start_id = 0
            non_basic_coeff_id = 0
            max_val = 0
            for i in range(self.n):

                if (self.C_j[self.non_basic_var[i]] - np.dot(self.y, np.transpose(self.A_ij[:, self.non_basic_var[i]]))) > 0:
                    next_loop = True

                if (self.C_j[self.non_basic_var[i]] - np.dot(self.y, np.transpose(self.A_ij[:, self.non_basic_var[i]]))) > max_val:
                    max_val = self.C_j[self.non_basic_var[i]] - np.dot(self.y, np.transpose(self.A_ij[:, self.non_basic_var[i]]))
                    start_id = self.non_basic_var[i]
                    non_basic_coeff_id = i
            if next_loop == False:
                break

            print(f'\n\nMax_j Cj - Zj = {max_val}')
            print(f'Entering variable: x_{start_id+1}')
            iter_number+=1
            if iter_number >= self.MAX_itr:
                print("\n\nCaught into infinite loop problem is unbounded\n")
                break
        pass

    def printsolution(self):
        print('\n\n||   Final solution   || ')
        Z_val = 0
        for i in range(self.m):
            print(f'x_{int(self.basic_var[i])+1} = {self.X_B[i]}')
            Z_val += self.C_j[int(self.basic_var[i])] * self.X_B[i]
        print(f'Z = {Z_val}')

    def run(self, file_name):

        if type(file_name) == type(None):
            self.input()
        else:
            self.file_input(file_name)

        start_id, nb_coeff_id = self.initilize()

        self.printtable()

        self.compute(start_id, nb_coeff_id)

        self.printsolution()
        pass


RS = Revised_Simplex()

RS.run('test2.txt') # uncomment for problem A 

# RS.run('B.txt') # uncomment for problem B

# RS.run('C.txt') # uncomment for problem C

# RS.run('D.txt') # uncomment for problem D

# RS.run('E.txt') # uncomment for problem E

# to run with manual input uncomment the line below and comment the line above

# RS.run()




