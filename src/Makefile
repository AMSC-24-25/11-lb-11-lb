# Variabili
CXX = g++                       # Compilatore
CXXFLAGS = -std=c++17 -O2 -Wall # Flag di compilazione
SRC_DIR = .                      # Directory dei sorgenti (cartella corrente)
TARGET = lbm_sim                 # Nome dell'eseguibile finale

# Trova i file sorgente .cpp
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)  # Trova tutti i file .cpp nella directory
OBJ_FILES = $(SRC_FILES:.cpp=.o)          # Crea i file .o corrispondenti ai file .cpp

# Regola principale
all: $(TARGET)

# Regola per compilare i file oggetto
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@  # Compila il file .cpp in un file .o

# Regola per creare l'eseguibile
$(TARGET): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(OBJ_FILES) -o $(TARGET)  # Linka i file .o per creare l'eseguibile

# Regola per pulire
clean:
	rm -f $(OBJ_FILES) $(TARGET)  # Elimina i file oggetto e l'eseguibile

# Targets fittizi
.PHONY: all clean

