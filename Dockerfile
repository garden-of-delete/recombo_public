FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    cmake \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /recombo
COPY . .

RUN cmake -B build && cmake --build build
RUN ./build/src/bin/unitTest
