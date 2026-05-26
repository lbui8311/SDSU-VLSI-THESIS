`timescale 1ns / 1ps

module priority_encoder #(parameter N = 8)
                         (input wire [N-1:0] in,
                          output reg [$clog2(N)-1:0] code, // Binary index of the selected input
                          output reg valid); // Asserted when any input bit is high
                          
    integer i; // For-loop incrementing index
    
    always @* begin
        if (in == 0) begin // Set both ouputs to zero if input is flat zero (all bits are zeros)
            valid = 0;
            code = 0;
        end else begin
            valid = 1; // If input is not flat zero, valid flag is set high
            code = 0;  // Default output if output is not set in for-loop
            begin : loop
                for (i = N-1; i >= 0; i = i - 1) begin
                    if (in[i] == 1'b1) begin
                        code = i;
                        disable loop; // break out of loop
                    end
                end
            end
        end
    end              
endmodule
