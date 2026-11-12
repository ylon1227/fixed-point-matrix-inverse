/* Minimal bare-metal startup for STM32F4 (no HAL, no CMSIS). License: MIT */
#include <stdint.h>

extern uint32_t _estack, _sidata, _sdata, _edata, _sbss, _ebss;
int main(void);

void Reset_Handler(void)
{
    uint32_t *src = &_sidata, *dst = &_sdata;
    while (dst < &_edata) *dst++ = *src++;
    for (dst = &_sbss; dst < &_ebss; ) *dst++ = 0;
    /* enable the FPU (CP10/CP11 full access) for the hardware-float reference */
    *(volatile uint32_t *)0xE000ED88 |= (0xFu << 20);
    __asm volatile ("dsb\n isb");
    main();
    for (;;) { }
}

void Default_Handler(void) { for (;;) { } }

__attribute__((section(".isr_vector"), used))
const void *const g_vectors[16] = {
    &_estack, (void *)Reset_Handler,
    (void *)Default_Handler, (void *)Default_Handler, (void *)Default_Handler,
    (void *)Default_Handler, (void *)Default_Handler, 0, 0, 0, 0,
    (void *)Default_Handler, (void *)Default_Handler, 0,
    (void *)Default_Handler, (void *)Default_Handler,
};
