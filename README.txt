./ECGDetectionSoftware - skrypt uruchamiaj¹cy aplikacjê
./ECGAnalyzer - plik zawieraj¹cy klasê ECGAnalyzer s³u¿¹c¹ do analizy sygna³u
./signalutils - plik zawieraj¹cy zestaw przydatnych funkcji
./Spectrum - plik zawieraj¹cy klasê przedstawiaj¹c¹ widmo amplitudowe sygna³u
./test - skrypt generuj¹cy szum do testów algorytmy
./test/ - katalog z zasobami do testów algorytmu
./test/test1/ - katalog zawieraj¹cy sygna³y EKG pacjentów wykorzystane w pierwszym teœcie
./test/test2/ - katalog zawieraj¹cy sygna³y EKG z szumem wykorzystane w drugim teœcie
./sampledata/ - katalog zawieraj¹cy sygna³y EKG pacjentów z bazy PhysioNet
./requirements - plik zawieraj¹cy listê potrzebnych pakietów do uruchomienia aplikacji

Instalacja Pythona i zale¿noœci:
- Instalacja Python3.6 [1] wraz z pip [2] (powinien byæ dostarczany z Pythonem3.6)
- Instalacja pakietów
$ pip3.6 install -r ./requirements.txt

Uruchamianie aplikacji:
$ python3.6 ./ECGDetectionSoftware

[1] https://www.python.org/downloads/
[2] https://pypi.python.org/pypi/pip