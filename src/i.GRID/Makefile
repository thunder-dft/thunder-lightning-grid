GRID = M_assemble_2c_grid.o M_assemble_usr_grid.o M_fft_grid.o M_project_grid.o M_density_matrix_grid.o

all: libthd_GRID.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

libthd_GRID.a: $(GRID)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
