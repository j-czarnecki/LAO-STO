#######################################
########### Dependencies ##############
#######################################

$(OBJ_DIR)/input_output/src/reader.o: $(OBJ_DIR)/physical/src/parameters.o \
				  $(OBJ_DIR)/input_output/src/logger.o \
				  $(OBJ_DIR)/types/types.o \
				  $(OBJ_DIR)/physical/src/utilities.o

$(OBJ_DIR)/input_output/src/logger.o:

$(OBJ_DIR)/input_output/src/macros_def.o:

$(OBJ_DIR)/input_output/src/writers.o: $(OBJ_DIR)/physical/src/parameters.o \
				   $(OBJ_DIR)/input_output/src/reader.o \
				   $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/physical/src/parameters.o:

$(OBJ_DIR)/types/types.o: $(OBJ_DIR)/physical/src/parameters.o

$(OBJ_DIR)/physical/src/utilities.o: $(OBJ_DIR)/physical/src/parameters.o \
				 $(OBJ_DIR)/types/types.o
