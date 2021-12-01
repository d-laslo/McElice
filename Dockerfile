FROM python:3.9-slim
COPY . .
RUN pip install numpy
CMD ["python3","McEliece.py"]
