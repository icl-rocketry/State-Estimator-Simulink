function connect_to_ricardo_udp()

    port = 5005;   % <-- make sure this matches backend
    disp("🔌 Setting up Python UDP client...");

    py.importlib.import_module('socket');

    sock = py.socket.socket(py.socket.AF_INET, py.socket.SOCK_DGRAM);
    sock.bind(py.tuple({'0.0.0.0', int32(port)}));
    sock.settimeout(1);

    disp("✅ Listening for UDP telemetry on port " + string(port));
    disp("📡 Waiting for packets...");

    while true
        try
            % recvfrom returns a Python tuple → MATLAB sees a cell array
            result = sock.recvfrom(int32(4096));
            msg  = result{1};
            addr = result{2};

            raw = char(msg);
            raw = regexprep(raw,'[^\d\.\-\s]',''); % remove non-numeric
            data = str2num(raw); %#ok<ST2NM>

            disp("📨 UDP Packet:");
            disp(data);

        catch ME
            if contains(string(ME.message),"timed out")
                disp("⏳ No UDP packet received...");
            else
                rethrow(ME);
            end
        end
    end
end
