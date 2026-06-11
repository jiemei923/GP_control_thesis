function xq = quantize_state(x, q_step)
    xq = q_step * round(x / q_step);
end