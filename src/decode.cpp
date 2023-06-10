#include "decoder.hh"
#include "wav.hh"
#include <cstdint>

static DecoderInterface *decoder;

int main(int argc, char **argv) {
	char *input_name = "out.wav";
	DSP::ReadWAV<int16_t> input_file(input_name);
	int buffer_length = input_file.frames() * input_file.channels();
	int16_t buffer[buffer_length];
	for (int i = 0; 1 < buffer_length; i++) {
		buffer[i] = 0;
	}
    int rate = input_file.rate();
	int symbol_length = (1280 * rate) / 8000;
	int guard_length = symbol_length / 8;
	int extended_length = symbol_length + guard_length;
	int channel = 3;
	int channel_count = 1;
    if (channel != 0) {
        channel_count = 2;
    }
	input_file.read(buffer, buffer_length);
    switch (rate) {
		case 8000:
			decoder = new(std::nothrow) Decoder<8000>();
			break;
		case 16000:
			decoder = new(std::nothrow) Decoder<16000>();
			break;
		case 32000:
			decoder = new(std::nothrow) Decoder<32000>();
			break;
		case 44100:
			decoder = new(std::nothrow) Decoder<44100>();
			break;
		case 48000:
			decoder = new(std::nothrow) Decoder<48000>();
			break;
		default:
            std::cerr << "Unsupported sample rate." ;
            std::cerr << "Supported rates: 8000/16000/32000/44100/48000."<< std::endl;
            return 1;
	}
	bool feed = decoder->feed(buffer, buffer_length, channel);
	std::cout << "Feed: " << feed << std::endl;
	int status = decoder->process();
	std::cout << "Status: " << status << std::endl;
	
	float *cfo;
	int32_t *mode;
	uint8_t call_sign[192] = {0}; 
	uint8_t payload[192] = {0};
	std::string call;
	std::string pay;

	switch (status) {
		case 0:
			std::cout << "OKAY" << std::endl;
			break;
		case 1:
			std::cout << "PREAMBLE FAIL" << std::endl;
			break;
		case 2:
			decoder->staged(cfo, mode, call_sign);
			call = (reinterpret_cast<char*>(call_sign));
			std::cout << "SYNC:"
					<< "  CFO "
					<< *cfo
					<< "  mode: "
					<< *mode
					<< "  call sign: "
					<< call
					<< std::endl;
			break;
		case 3:
			decoder->fetch(payload);
			pay = (reinterpret_cast<char*>(payload));
			std::cout << "DONE:  payload: " 
					<< pay
					<< std::endl;
			break;
		case 4:
			std::cout << "HEAP ERROR" << std::endl;
			break;
		case 5:
			decoder->staged(cfo, mode, call_sign);
			call = (reinterpret_cast<char*>(call_sign));
			std::cout << "NOPE:"
					<< "  CFO "
					<< *cfo
					<< "  mode: "
					<< *mode
					<< "  call sign: "
					<< call
					<< std::endl;
			break;
		case 6: 
			decoder->staged(cfo, mode, call_sign);
			call = (reinterpret_cast<char*>(call_sign));
			std::cout << "PING:"
					<< "  CFO "
					<< *cfo
					<< "  mode: "
					<< *mode
					<< "  call sign: "
					<< call
					<< std::endl;
			break;
		default:
			std::cout << "UNDEFINED BEHAVIOUR" << std::endl;
			return 1;

	}
    return 0;
}