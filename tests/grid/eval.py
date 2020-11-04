import numpy as np
from time import time

# Custom library
import _gridprocessing

def run():
    print("\n[==> Creating grids")

    print("\n> Integer grid")
    start_int = time()
    custom_grid = _gridprocessing.igrid(3 * 4 * 4).reshape(3, 4, 4)
    end_int = time()
    print(custom_grid)
    # Calculate time
    time_int = end_int - start_int

    print("\n> Float grid")
    start_float = time()
    custom_grid = _gridprocessing.fgrid(3 * 4 * 4).reshape(3, 4, 4)
    end_float = time()
    print(custom_grid)
    # Calculate time
    time_float = end_float - start_float

    print("\n> Double grid")
    start_double = time()
    custom_grid = _gridprocessing.dgrid(3 * 4 * 4).reshape(3, 4, 4)
    print(custom_grid)
    end_double = time()
    # Calculate time
    time_double = end_double - start_double

    print("\n> Char grid")
    start_char = time()
    custom_grid = _gridprocessing.cgrid(3 * 4 * 4).reshape(3, 4, 4)
    print(custom_grid)
    end_char = time()
    # Calculate time
    time_char = end_char - start_char

    print("\n> Numpy integer grid")
    start_np_int = time()
    numpy_grid = np.ones( [3, 4, 4] , dtype=np.intc)
    print(numpy_grid)
    end_np_int = time()
    # Calculate time
    time_np_int = end_np_int - start_np_int

    print("\n> Numpy float grid")
    start_np_float = time()
    numpy_grid = np.zeros( [3, 4, 4] , dtype=np.single)
    print(numpy_grid)
    end_np_float = time()
    # Calculate time
    time_np_float = end_np_float - start_np_float

    print("\n> Numpy double grid")
    start_np_double = time()
    numpy_grid = np.zeros( [3, 4, 4] , dtype=np.double)
    print(numpy_grid)
    end_np_double = time()
    # Calculate time
    time_np_double = end_np_double - start_np_double

    print("\n> Numpy char grid")
    start_np_char = time()
    numpy_grid = np.char.chararray( [3, 4, 4] )
    print(numpy_grid)
    end_np_char = time()
    # Calculate time
    time_np_char = end_np_char - start_np_char

    print("\n\t\tTimes")
    print(f"Integer Grid: \t\t{time_int:0.6f}")
    print(f"Float Grid: \t\t{time_float:0.6f}")
    print(f"Double Grid: \t\t{time_double:0.6f}")
    print(f"Char Grid: \t\t{time_char:0.6f}")
    print(f"Numpy integer Grid: \t{time_np_int:0.6f}")
    print(f"Numpy float Grid: \t{time_np_float:0.6f}")
    print(f"Numpy double Grid: \t{time_np_double:0.6f}")
    print(f"Numpy char Grid: \t{time_np_char:0.6f}")

