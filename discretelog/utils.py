from subprocess import Popen, PIPE, STDOUT
from multiprocessing import Process, Queue, cpu_count
from typing import Callable, Iterable, Optional, Any
from tqdm import trange, tqdm


def mrange(
    a: int, b: Optional[int] = None, c: int = 1, DEBUG: bool = False
) -> Iterable[int]:
    if DEBUG:
        if b is None:
            return trange(a, leave=False)
        else:
            return trange(a, b, c, leave=False)
    else:
        if b is None:
            return range(a)
        else:
            return range(a, b, c)


def mexec(command: str | list[str], DEBUG: bool = False) -> str:
    if DEBUG:
        print(f"Executing {command}")
    if isinstance(command, str):
        command = command.split()
    captured_output = []
    with Popen(command, stdout=PIPE, stderr=STDOUT, text=True) as process:
        while True:
            if process.stdout is not None:
                output_line = process.stdout.readline()
                if not output_line and process.poll() is not None:
                    break
                if DEBUG:
                    print(output_line, end="")
                captured_output.append(output_line)
            else:
                break
    return "".join(captured_output).strip()


def parallel_for(
    f: Callable[[Any], Any],
    params: tuple[int, ...],
    prange: range,
    cores: Optional[int] = None,
    DEBUG: bool = False,
) -> dict[int, int]:
    n = len(prange)
    workers = cpu_count() if cores is None else cores

    pb = tqdm(total=n) if DEBUG else NoOpTqdm()
    pb.refresh()

    it = iter(prange)
    procs: list[Process] = []
    res = {}
    started, ended = 0, 0
    queue: Queue[tuple[int, int]] = Queue()
    while n > ended:
        if n > started and len(procs) - ended < workers:
            started += 1
            j = next(it)
            p = Process(target=f, args=(queue, j, *params))
            p.start()
            procs += [p]
        while not queue.empty():
            j, resj = queue.get()
            res[j] = resj
            ended += 1

            pb.update(1)

    pb.close()

    for p in procs:
        p.join()
        p.close()

    return res


def parallel_for_balanced(
    f: Callable[..., Any],
    params: tuple[Any, ...],
    prange: range,
    cores: Optional[int] = None,
    DEBUG: bool = False,
) -> dict[int, Any]:
    n = len(prange)
    workers = cpu_count() if cores is None else cores
    pb = tqdm(total=n) if DEBUG else NoOpTqdm()
    pb.refresh()
    its: list[list[int]] = [[] for _ in range(workers)]
    for i, it in enumerate(prange):
        its[i % workers].append(it)
    res = {}
    procs: list[Process] = []
    queue: Queue[tuple[int, int]] = Queue()
    for i in range(workers):
        p = Process(target=f, args=(queue, its[i], *params))
        p.start()
        procs += [p]
    ended = 0
    while ended < workers:
        while not queue.empty():
            c = queue.get()
            if c is not None:
                j, resj = c
                res[j] = resj
                if DEBUG:
                    pb.update(1)
            else:
                ended += 1
    for p in procs:
        p.join()
        p.close()
    return res


class NoOpTqdm:
    def update(self, _: int) -> None:
        pass

    def close(self) -> None:
        pass

    def refresh(self) -> None:
        pass
