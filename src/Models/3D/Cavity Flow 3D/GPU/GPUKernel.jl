using CUDA

function add_broadcast!(y, x)
    CUDA.@sync y .+= x
    return
end

function substract_broadcast!(y, x)
    CUDA.@sync y .-= x
    return
end

function product_broadcast!(y, x)
    CUDA.@sync y .-= x
    return
end