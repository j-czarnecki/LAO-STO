#######################################
########### Dependencies ##############
#######################################

$(OBJ_DIR)/physical/src/utilities.o: $(OBJ_DIR)/physical/src/parameters.o \
						 $(OBJ_DIR)/input_output/reader.o

$(OBJ_DIR)/physical/src/parameters.o:

$(OBJ_DIR)/input_output/types.o: $(OBJ_DIR)/physical/src/parameters.o

$(OBJ_DIR)/physical/src/hamiltonians.o: $(OBJ_DIR)/physical/src/utilities.o \
								$(OBJ_DIR)/physical/src/parameters.o \
								$(OBJ_DIR)/input_output/reader.o \
								$(OBJ_DIR)/types/types.o

$(OBJ_DIR)/input_output/reader.o: $(OBJ_DIR)/physical/src/parameters.o \
					$(OBJ_DIR)/input_output/logger.o \
					$(OBJ_DIR)/types/types.o

$(OBJ_DIR)/input_output/logger.o:

$(OBJ_DIR)/input_output/macros_def.o:

$(OBJ_DIR)/input_output/writers.o: $(OBJ_DIR)/physical/src/parameters.o \
					$(OBJ_DIR)/input_output/reader.o \
					$(OBJ_DIR)/types/types.o