module LFSR8_Flipped(
    input wire clk,            // Clock input
    input wire reset,          // Reset signal
    output reg [7:0] LFSR_out, // LFSR normal output
    output reg [7:0] LFSR_flip // LFSR flipped output
);

    reg [7:0] LFSR; // Register to hold the LFSR state

    always @(posedge clk or posedge reset) begin
        if (reset) begin
            // Initialize LFSR to 8'b10000000, cannot be zero
            LFSR <= 8'b10000000;
        end else begin
            // Update LFSR with feedback logic
            LFSR <= {xor(xor(LFSR[7], LFSR[5]), xor(LFSR[4], LFSR[3])), LFSR[7:1]};
        end
    end

    // Assign outputs
    always @(posedge clk or posedge reset) begin
        LFSR_out <= LFSR;       // Normal LFSR output
        LFSR_flip <= {LFSR[0], LFSR[1], LFSR[2], LFSR[3], LFSR[4], LFSR[5], LFSR[6], LFSR[7]}; // Flipped output
    end
endmodule
