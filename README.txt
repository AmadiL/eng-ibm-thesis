./ECGDetectionSoftware - skrypt uruchamiaj�cy aplikacj�
./ECGAnalyzer - plik zawieraj�cy klas� ECGAnalyzer s�u��c� do analizy sygna�u
./signalutils - plik zawieraj�cy zestaw przydatnych funkcji
./Spectrum - plik zawieraj�cy klas� przedstawiaj�c� widmo amplitudowe sygna�u
./test - skrypt generuj�cy szum do test�w algorytmy
./test/ - katalog z zasobami do test�w algorytmu
./test/test1/ - katalog zawieraj�cy sygna�y EKG pacjent�w wykorzystane w pierwszym te�cie
./test/test2/ - katalog zawieraj�cy sygna�y EKG z szumem wykorzystane w drugim te�cie
./sampledata/ - katalog zawieraj�cy sygna�y EKG pacjent�w z bazy PhysioNet
./requirements - plik zawieraj�cy list� potrzebnych pakiet�w do uruchomienia aplikacji

Instalacja Pythona i zale�no�ci:
- Instalacja Python3.6 [1] wraz z pip [2] (powinien by� dostarczany z Pythonem3.6)
- Instalacja pakiet�w
$ pip3.6 install -r ./requirements.txt

Uruchamianie aplikacji:
$ python3.6 ./ECGDetectionSoftware

[1] https://www.python.org/downloads/
[2] https://pypi.python.org/pypi/pip