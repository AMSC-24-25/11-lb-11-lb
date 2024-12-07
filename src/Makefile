# Variabili
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall
SRC_DIR = .
TARGET = lbm_sim

# Trova tutti i file sorgente .cpp e genera i corrispondenti file oggetto .o
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(SRC_FILES:.cpp=.o)

# Regola principale
all: $(TARGET)

# Regola per creare l'eseguibile
$(TARGET): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(OBJ_FILES) -o $@

# Regola per compilare ogni file oggetto
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Regola per pulire i file generati
clean:
	rm -f $(OBJ_FILES) $(TARGET)

# Target fittizi
.PHONY: all clean