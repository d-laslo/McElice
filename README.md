# McEliece

To run in docker use this commands

```
docker build -t mceliece .
docker run mceliece
```

To run the raw test

```
python3 testMcEliece.py
```

To run the profilling. The profilling result will be written to file profile.txt

```
python3 -m cProfile testMcEliece.py > profile.txt
```
