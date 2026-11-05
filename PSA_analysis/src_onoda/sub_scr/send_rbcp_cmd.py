import socket
import struct


# RBCP定数
RBCP_VER = 0xFF
RBCP_CMD_WR = 0x80
RBCP_CMD_RD = 0xC0
RBCP_ACK = 0x08

RBCP_HEADER_SIZE = 8
UDP_BUF_SIZE = 2048

DEFAULT_IP = "192.168.100.128"
DEFAULT_PORT = 4660


def format_hex(data):
    return " ".join(f"{byte:02X}" for byte in data)


def parse_rbcp_response(response):
    """RBCPレスポンスをヘッダとデータに分解する。"""
    if response is None:
        return None

    if len(response) < RBCP_HEADER_SIZE:
        raise ValueError(
            f"RBCP response is too short: {len(response)} bytes"
        )

    version, command, packet_id, length, address = struct.unpack(
        "!BBBBI",
        response[:RBCP_HEADER_SIZE]
    )

    data = response[RBCP_HEADER_SIZE:]

    if len(data) != length:
        raise ValueError(
            f"Invalid data length: header={length}, actual={len(data)}"
        )

    return {
        "version": version,
        "command": command,
        "packet_id": packet_id,
        "length": length,
        "address": address,
        "data": data,
    }


def print_rbcp_response(response):
    """RBCPレスポンスを読みやすく表示する。"""
    parsed = parse_rbcp_response(response)

    if parsed is None:
        print("No response")
        return

    print(
        f"RBCP response: "
        f"version=0x{parsed['version']:02X}, "
        f"command=0x{parsed['command']:02X}, "
        f"id=0x{parsed['packet_id']:02X}, "
        f"length={parsed['length']}, "
        f"address=0x{parsed['address']:08X}"
    )

    if parsed["data"]:
        print(f"Data: {format_hex(parsed['data'])}")

def send_rbcp_command(
    ip,
    port,
    command,
    address,
    length,
    data=b"",
    verbose=False
):
    if not 0 <= length <= 0xFF:
        raise ValueError("RBCP length must be between 0 and 255")

    if command == RBCP_CMD_WR and len(data) != length:
        raise ValueError(
            f"Write data length mismatch: length={length}, data={len(data)}"
        )

    if not 0 <= address <= 0xFFFFFFFF:
        raise ValueError(
            "RBCP address must be between 0x00000000 and 0xFFFFFFFF"
        )

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.settimeout(1.0)

    header = struct.pack(
        "!BBBBI",
        RBCP_VER,
        command,
        0x01,
        length,
        address
    )

    packet = header + data

    try:
        if verbose:
            print(f"send to  : {ip}:{port}")
            print(f"command  : {format_hex(packet)}")

        sock.sendto(packet, (ip, port))

        response, sender = sock.recvfrom(UDP_BUF_SIZE)

        if verbose:
            print(f"recv from: {sender[0]}:{sender[1]}")
            print(f"response : {format_hex(response)}")

        return response

    except socket.timeout:
        if verbose:
            print("response : Timeout")

        return None

    finally:
        sock.close()


def parse_and_execute(
    cmdline,
    ip=DEFAULT_IP,
    port=DEFAULT_PORT,
    verbose=False
):
    tokens = cmdline.strip().split()

    if not tokens:
        return None

    cmd = tokens[0].lower()

    try:
        if cmd == "wrb":
            if len(tokens) != 3:
                raise ValueError("Usage: wrb [addr] [val]")

            addr = int(tokens[1], 0)
            value = int(tokens[2], 0)

            if not 0 <= value <= 0xFF:
                raise ValueError("wrb value must be between 0x00 and 0xFF")

            data = struct.pack("!B", value)

            return send_rbcp_command(
                ip,
                port,
                RBCP_CMD_WR,
                addr,
                1,
                data,
                verbose=verbose
            )

        elif cmd == "wrs":
            if len(tokens) != 3:
                raise ValueError("Usage: wrs [addr] [val]")

            addr = int(tokens[1], 0)
            value = int(tokens[2], 0)

            if not 0 <= value <= 0xFFFF:
                raise ValueError("wrs value must be between 0x0000 and 0xFFFF")

            data = struct.pack("!H", value)

            return send_rbcp_command(
                ip,
                port,
                RBCP_CMD_WR,
                addr,
                2,
                data,
                verbose=verbose
            )

        elif cmd == "wrw":
            if len(tokens) != 3:
                raise ValueError("Usage: wrw [addr] [val]")

            addr = int(tokens[1], 0)
            value = int(tokens[2], 0)

            if not 0 <= value <= 0xFFFFFFFF:
                raise ValueError(
                    "wrw value must be between 0x00000000 and 0xFFFFFFFF"
                )

            data = struct.pack("!I", value)

            return send_rbcp_command(
                ip,
                port,
                RBCP_CMD_WR,
                addr,
                4,
                data,
                verbose=verbose
            )

        elif cmd == "rd":
            if len(tokens) != 3:
                raise ValueError("Usage: rd [addr] [len]")

            addr = int(tokens[1], 0)
            length = int(tokens[2], 0)

            if not 1 <= length <= 0xFF:
                raise ValueError("Read length must be between 1 and 255")

            return send_rbcp_command(
                ip,
                port,
                RBCP_CMD_RD,
                addr,
                length,
                b"",
                verbose=verbose
            )

        elif cmd == "load":
            if len(tokens) != 2:
                raise ValueError("Usage: load [file]")

            filename = tokens[1]
            response = None

            with open(filename, encoding="utf-8") as file:
                for line_number, line in enumerate(file, start=1):
                    line = line.strip()

                    # 空行とコメント行を無視
                    if not line or line.startswith("#"):
                        continue

                    try:
                        response = parse_and_execute(
                            line,
                            ip,
                            port,
                            verbose=verbose
                        )
                    except (ValueError, OSError) as exc:
                        print(
                            f"{filename}:{line_number}: {exc}"
                        )

            return response

        elif cmd == "help":
            print("wrb [addr] [val]")
            print("wrs [addr] [val]")
            print("wrw [addr] [val]")
            print("rd  [addr] [len]")
            print("load [file]")
            print("quit")
            return None

        elif cmd in ("quit", "exit"):
            raise SystemExit

        else:
            print(f"Unknown command: {cmd}")
            return None

    except (IndexError, ValueError) as exc:
        print(f"Command error: {exc}")
        return None


def main():
    import sys

    if len(sys.argv) == 3:
        ip = sys.argv[1]

        try:
            port = int(sys.argv[2], 0)
        except ValueError:
            print(f"Invalid port: {sys.argv[2]}")
            return
    else:
        ip = DEFAULT_IP
        port = DEFAULT_PORT

    while True:
        try:
            cmdline = input("SiTCP-RBCP$ ")
            parse_and_execute(
                cmdline,
                ip,
                port,
                verbose=True
            )

        except (EOFError, KeyboardInterrupt):
            print("\nExiting.")
            break

        except SystemExit:
            print("Exiting.")
            break


if __name__ == "__main__":
    main()