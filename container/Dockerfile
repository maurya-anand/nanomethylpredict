FROM julia:1.12.1

WORKDIR /app

RUN apt-get update && apt-get install -y \
    curl unzip wget gzip \
    samtools bedtools bcftools tabix && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV JULIA_PKG_SERVER=https://mirrors.pku.edu.cn/julia

RUN curl -L "https://gitee.com/zhixingfeng/NanoFreeLunch.jl/archive/refs/tags/0.28.0.zip" -o nanofreelunch.zip && \
    unzip nanofreelunch.zip && \
    rm nanofreelunch.zip

WORKDIR /app/NanoFreeLunch.jl-0.28.0

RUN sed -i '761s/@main/Comonicon.@main/' src/NanoFreeLunch.jl

RUN julia setup.jl

ENV PATH="/app/NanoFreeLunch.jl-0.28.0/build/bin:$PATH"
